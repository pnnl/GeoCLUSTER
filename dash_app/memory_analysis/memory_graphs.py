import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import re

def parse_profile_output(file_path):
    data = []
    current_file = None
    function_start = True
    sequence_number = 0
    
    # Try different encodings
    # encodings = ['utf-8-sig', 'utf-8', 'latin1', 'cp1252']
    encoding = 'latin1'
    try:
        with open(file_path, 'r', encoding=encoding) as f:
            for line in f:
                # Check for function name and reset function_start flag
                if 'Filename:' in line:
                    function_match = re.search(r'Filename: .*\\([^\\]+\.py)', line)
                    if function_match:
                        current_file = function_match.group(1)
                        function_start = True
                    continue
                
                # Parse memory usage lines
                mem_match = re.match(r'\s*(\d+)\s+([\d.]+)\sMiB\s+([-\d.]+)\sMiB\s+(\d+)\s+(.*)', line)
                if mem_match:
                    line_num, mem_usage, increment, occurrences, line_content = mem_match.groups()
                    sequence_number += 1
                    
                    # Skip the first memory increase for each function
                    if function_start and float(increment) > 0:
                        function_start = False
                        continue
                    
                    data.append({
                        'File': current_file,
                        'Line': int(line_num),
                        'Sequence': sequence_number,
                        'Memory_Usage': float(mem_usage),
                        'Increment': float(increment),
                        'Occurrences': int(occurrences),
                        'Content': line_content.strip(),
                        'Label': f"{current_file}:{line_num}"
                    })
    except Exception as e:
        print(f"Error reading file with encoding {encoding}: {str(e)}")
        return pd.DataFrame(columns=['File', 'Line', 'Sequence', 'Memory_Usage', 
                                  'Increment', 'Occurrences', 'Content', 'Label'])

    if not data:
        print("No data was parsed from the file")
        return pd.DataFrame(columns=['File', 'Line', 'Sequence', 'Memory_Usage', 
                                  'Increment', 'Occurrences', 'Content', 'Label'])
    
    return pd.DataFrame(data)

def plot_individual_files(file_path):
    df = parse_profile_output(file_path)
    
    if df.empty:
        print("No data to plot")
        return
    
    unique_files = df['File'].unique()
    
    if len(unique_files) == 0:
        print("No files found in the data")
        return
    
    # Create a figure for each file
    for file in unique_files:
        file_data = df[df['File'] == file]
        
        # Create a figure with 2 subplots side by side
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        fig.suptitle(f'Memory Profile for {file}', fontsize=16)
        
        # Plot 1: Memory Usage Over Time
        ax1.plot(file_data['Sequence'], file_data['Memory_Usage'], 
                marker='o', markersize=4, linewidth=2)
        ax1.set_xlabel('Sequential Step')
        ax1.set_ylabel('Memory Usage (MiB)')
        ax1.set_title('Memory Usage Over Time')
        ax1.grid(True)
        
        # Add annotations for significant memory changes
        significant_changes = file_data[abs(file_data['Increment']) > 0.1]
        for _, row in significant_changes.iterrows():
            ax1.annotate(f"+{row['Increment']:.1f} MiB\nLine {row['Line']}",
                        (row['Sequence'], row['Memory_Usage']),
                        xytext=(10, 10), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                        arrowprops=dict(arrowstyle='->'))
        
        # Plot 2: Memory Increments
        significant_increments = file_data[file_data['Increment'] > 0.1].copy()
        if not significant_increments.empty:
            bars = ax2.bar(significant_increments['Sequence'], 
                         significant_increments['Increment'])
            
            # Add value labels on top of bars
            for bar in bars:
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height,
                        f'{height:.1f}',
                        ha='center', va='bottom')
            
            # Add content annotations for significant increments
            for _, row in significant_increments.iterrows():
                ax2.annotate(f"Line {row['Line']}: {row['Content'][:30]}...",
                            (row['Sequence'], row['Increment']),
                            xytext=(10, 10), textcoords='offset points',
                            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.5),
                            arrowprops=dict(arrowstyle='->'))
                
            ax2.set_xlabel('Sequential Step')
            ax2.set_ylabel('Memory Increment (MiB)')
            ax2.set_title('Significant Memory Increments (>0.1 MiB)')
        else:
            ax2.text(0.5, 0.5, 'No significant memory increments found',
                    ha='center', va='center', transform=ax2.transAxes)
        
        plt.tight_layout()
        plt.show()
        
        # Print detailed summary for this file
        print(f"\nDetailed Summary for {file}")
        print("-" * 80)
        print(f"Total memory increase: {file_data['Increment'].sum():.2f} MiB")
        print(f"Maximum memory usage: {file_data['Memory_Usage'].max():.2f} MiB")
        print(f"Number of significant increases (>0.1 MiB): {len(significant_increments)}")
        
        if not significant_increments.empty:
            print("\nTop memory-consuming lines:")
            top_lines = significant_increments.nlargest(5, 'Increment')
            for _, row in top_lines.iterrows():
                print(f"\nStep {row['Sequence']} (Line {row['Line']}): {row['Increment']:.2f} MiB")
                print(f"Content: {row['Content'][:100]}...")
        print("-" * 80)

def plot_memory_usage(file_path):
    df = parse_profile_output(file_path)
    
    if df.empty:
        print("No data to plot")
        return
    
    plt.figure(figsize=(15, 6))
    
    # Memory usage over sequential steps with file information
    plt.subplot(1, 2, 1)
    for file in df['File'].unique():
        file_data = df[df['File'] == file]
        plt.plot(file_data['Sequence'], file_data['Memory_Usage'], 
                label=file, marker='o', markersize=4)
    
    plt.xlabel('Sequential Step')
    plt.ylabel('Memory Usage (MiB)')
    plt.title('Memory Usage Over Time\n(Excluding Initial Function Loads)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)

    # Memory increments
    plt.subplot(1, 2, 2)
    significant_increments = df[df['Increment'] > 0.1].copy()
    
    if not significant_increments.empty:
        colors = sns.color_palette("husl", n_colors=len(df['File'].unique()))
        file_color_dict = dict(zip(df['File'].unique(), colors))
        
        # Create labels that include sequence number
        significant_increments['SequenceLabel'] = significant_increments.apply(
            lambda row: f"{row['File']}\nStep {row['Sequence']}\nLine {row['Line']}", axis=1)
        
        bars = plt.bar(range(len(significant_increments)), 
                      significant_increments['Increment'],
                      color=[file_color_dict[file] for file in significant_increments['File']])
        
        plt.xticks(range(len(significant_increments)), 
                  significant_increments['SequenceLabel'],
                  rotation=45, ha='right')
        
        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.1f}',
                    ha='center', va='bottom')
        
        plt.xlabel('File and Step Information')
        plt.ylabel('Memory Increment (MiB)')
        plt.title('Significant Memory Increments\n(>0.1 MiB)')
        
        # Add file-based legend
        handles = [plt.Rectangle((0,0),1,1, color=color) for color in file_color_dict.values()]
        plt.legend(handles, file_color_dict.keys(), 
                  bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        plt.text(0.5, 0.5, 'No significant memory increments found',
                ha='center', va='center')
    
    plt.tight_layout()
    plt.show()

def plot_memory_detailed(file_path):
    df = parse_profile_output(file_path)
    
    if df.empty:
        print("No data to plot")
        return
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 15))
    
    # Cumulative memory usage by file
    for file in df['File'].unique():
        file_data = df[df['File'] == file]
        ax1.plot(file_data['Sequence'], file_data['Memory_Usage'], 
                label=file, marker='o', markersize=4)
    
    ax1.set_title('Cumulative Memory Usage by File\n(Excluding Initial Function Loads)')
    ax1.set_xlabel('Sequential Step')
    ax1.set_ylabel('Memory Usage (MiB)')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(True)
    
    # Memory increments with file information
    significant_increments = df[df['Increment'] > 0.1].copy()
    if not significant_increments.empty:
        colors = sns.color_palette("husl", n_colors=len(df['File'].unique()))
        file_color_dict = dict(zip(df['File'].unique(), colors))
        
        # Create labels that include sequence number
        significant_increments['SequenceLabel'] = significant_increments.apply(
            lambda row: f"{row['File']}\nStep {row['Sequence']}\nLine {row['Line']}", axis=1)
        
        bars = ax2.bar(range(len(significant_increments)), 
                      significant_increments['Increment'],
                      color=[file_color_dict[file] for file in significant_increments['File']])
        
        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.1f}',
                    ha='center', va='bottom')
        
        ax2.set_xticks(range(len(significant_increments)))
        ax2.set_xticklabels(significant_increments['SequenceLabel'], rotation=45, ha='right')
        ax2.set_title('Significant Memory Increments (>0.1 MiB)')
        ax2.set_xlabel('File and Step Information')
        ax2.set_ylabel('Memory Increment (MiB)')
        
        # Add file-based legend
        handles = [plt.Rectangle((0,0),1,1, color=color) for color in file_color_dict.values()]
        ax2.legend(handles, file_color_dict.keys(), 
                  bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax2.text(0.5, 0.5, 'No significant memory increments found',
                ha='center', va='center', transform=ax2.transAxes)
    
    # Top memory consumers
    if not significant_increments.empty:
        top_consumers = significant_increments.nlargest(10, 'Increment')
        
        # Create detailed labels for top consumers
        top_consumers['DetailLabel'] = top_consumers.apply(
            lambda row: f"{row['File']}\nStep {row['Sequence']}, Line {row['Line']}", axis=1)
        
        bars = ax3.barh(range(len(top_consumers)), 
                       top_consumers['Increment'],
                       color=[file_color_dict[file] for file in top_consumers['File']])
        
        ax3.set_title('Top 10 Memory-Consuming Steps')
        ax3.set_xlabel('Memory Increment (MiB)')
        ax3.set_yticks(range(len(top_consumers)))
        ax3.set_yticklabels(top_consumers['DetailLabel'])
        
        # Add content annotations
        for i, row in enumerate(top_consumers.itertuples()):
            ax3.text(row.Increment, i, f"  {row.Content[:50]}...", 
                    va='center', fontsize=8)
    else:
        ax3.text(0.5, 0.5, 'No significant memory increments found',
                ha='center', va='center', transform=ax3.transAxes)
    
    plt.tight_layout()
    plt.show()

    # Print detailed summary
    print("\nDetailed Memory Usage Summary")
    print("-" * 80)
    print(f"Total memory increase: {df['Increment'].sum():.2f} MiB")
    print(f"Maximum memory usage: {df['Memory_Usage'].max():.2f} MiB")
    print(f"Number of significant increases (>0.1 MiB): {len(significant_increments)}")
    
    print("\nBreakdown by File:")
    for file in df['File'].unique():
        file_data = df[df['File'] == file]
        file_increments = file_data[file_data['Increment'] > 0.1]
        print(f"\n{file}:")
        print(f"  Total increase: {file_data['Increment'].sum():.2f} MiB")
        print(f"  Significant increases: {len(file_increments)}")
        print(f"  Max memory: {file_data['Memory_Usage'].max():.2f} MiB")
    
    if not significant_increments.empty:
        print("\nTop 5 Memory-Consuming Steps:")
        print("-" * 80)
        top_5 = significant_increments.nlargest(5, 'Increment')
        for _, row in top_5.iterrows():
            print(f"\nFile: {row['File']}")
            print(f"Step {row['Sequence']} (Line {row['Line']}): {row['Increment']:.2f} MiB")
            print(f"Content: {row['Content'][:100]}...")
    print("-" * 80)

def print_memory_summary(file_path):
    df = parse_profile_output(file_path)
    significant_increments = df[df['Increment'] > 0.1]
    
    print("Memory Usage Summary:")
    print("-" * 80)
    
    # Summary by file
    print("\nMemory Usage by File:")
    for file in df['File'].unique():
        file_data = df[df['File'] == file]
        print(f"\n{file}:")
        print(f"Total memory increase: {file_data['Increment'].sum():.2f} MiB")
        print(f"Significant increases (>0.1 MiB): {len(file_data[file_data['Increment'] > 0.1])}")
    
    print("\nTop 5 Memory-Consuming Lines:")
    print("-" * 80)
    
    top_5 = significant_increments.nlargest(5, 'Increment')
    for _, row in top_5.iterrows():
        print(f"File: {row['File']}")
        print(f"Line {row['Line']}: {row['Increment']:.2f} MiB")
        print(f"Content: {row['Content'][:100]}...")
        print("-" * 80)



def main():
    file_path = "./memory_profile_output.txt"
    
    try:
        # individual file plots
        print("\nGenerating individual file plots...")
        plot_individual_files(file_path)

        # overall memory plots
        plot_memory_usage(file_path)
        plot_memory_detailed(file_path)
        print_memory_summary(file_path)
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()
