#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import pandas as pd

def write_excelsheet(df_summary, df_subsurf_res_mass, df_subsurf_res_time, df_econ, geoCLUSTER_results_pathname, model=None, sbt_params=None):

    # -----------------------------------------------------------------------------------------------------------
    # Creates an Excel file comprised of all results.
    # -----------------------------------------------------------------------------------------------------------

    writer = pd.ExcelWriter(geoCLUSTER_results_pathname, engine='xlsxwriter') 

    df_summary = pd.DataFrame.from_dict(df_summary, orient='index').transpose()
    df_subsurf_res_mass = pd.DataFrame(df_subsurf_res_mass)
    df_subsurf_res_time = pd.DataFrame(df_subsurf_res_time)
    df_econ = pd.DataFrame.from_dict(df_econ, orient='index').transpose()

    # Remove "Select Interpolation" row from download (but keep it in on-screen table)
    if 'Variable Parameter' in df_summary.columns:
        df_summary = df_summary[df_summary['Variable Parameter'] != 'Select Interpolation'].reset_index(drop=True)

    df_summary.to_excel(writer, sheet_name='Summary', index=False, index_label=False, header=True)
    df_subsurf_res_mass.to_excel(writer, sheet_name='Subsurface Results - Mass Flow',index=False, index_label=False, header=True)
    df_subsurf_res_time.to_excel(writer, sheet_name='Subsurface Results - Time',index=False, index_label=False, header=True)
    df_econ.to_excel(writer, sheet_name='Economic Results',index=False, index_label=False, header=True)
    
	#################### Next Tab ##########################
    wrksheet = writer.sheets['Summary']
    wrksheet.set_column(0, 0, 30)
    wrksheet.set_column(1, 1, 20)
    wrksheet.set_column(2, 2, 35)
    wrksheet.set_column(3, 3, 10)
    wrksheet.set_column(4, 4, 40)
    wrksheet.set_column(5, 5, 10)

	#################### Next Tab ##########################
    wrksheet = writer.sheets['Subsurface Results - Mass Flow']
    wrksheet.set_column(0, 0, 20)
    wrksheet.set_column(1, 1, 30)
    wrksheet.set_column(2, 2, 45)
    wrksheet.set_column(3, 3, 30)
    wrksheet.set_column(4, 4, 45)

	#################### Next Tab ##########################
    wrksheet = writer.sheets['Subsurface Results - Time']
    wrksheet.set_column(0, 0, 10)
    wrksheet.set_column(1, 1, 20)
    wrksheet.set_column(2, 2, 30)
    wrksheet.set_column(3, 3, 30)
    wrksheet.set_column(4, 4, 20)
    wrksheet.set_column(5, 5, 30)
    wrksheet.set_column(6, 6, 30)

    #################### Next Tab ##########################
    wrksheet = writer.sheets['Economic Results']
    wrksheet.set_column(0, 0, 10)
    wrksheet.set_column(1, 1, 25)
    wrksheet.set_column(2, 2, 30)
    wrksheet.set_column(3, 3, 25)
    wrksheet.set_column(4, 4, 30)
    wrksheet.set_column(5, 5, 30)
    wrksheet.set_column(6, 6, 35)
    wrksheet.set_column(7, 7, 30)
    wrksheet.set_column(8, 8, 35)

    # Add SBT-specific parameters sheet if model is SBT
    if model in ["SBT V1.0", "SBT V2.0"] and sbt_params:
        _add_sbt_parameters_sheet(writer, model, sbt_params)

    writer.close()


def _add_sbt_parameters_sheet(writer, model, sbt_params):
    """
    Add a sheet with SBT-specific parameters organized as editable parameters and assumptions.
    """
    editable_params = []
    assumptions = []
    
    def format_val(key, val):
        if val == "-" or val is None:
            return "-"
        if key in ["Tsurf", "c_m", "rho_m", "Diameter1", "Diameter2", "pipe_roughness"]:
            if isinstance(val, (int, float)):
                return val
        elif key in ["mesh", "accuracy", "n_laterals", "HyperParam1", "HyperParam3", "HyperParam5"]:
            if isinstance(val, (int, float)):
                return val
        elif key == "lateral_multiplier":
            if isinstance(val, (int, float)):
                return val
        return val
    
    if sbt_params.get("Tsurf") not in [None, "-"]:
        editable_params.append(["Surface Temperature", format_val("Tsurf", sbt_params.get("Tsurf")), "°C"])
    if sbt_params.get("c_m") not in [None, "-"]:
        editable_params.append(["Rock Specific Heat Capacity", format_val("c_m", sbt_params.get("c_m")), "J/kg-°C"])
    if sbt_params.get("rho_m") not in [None, "-"]:
        editable_params.append(["Rock Density", format_val("rho_m", sbt_params.get("rho_m")), "kg/m³"])
    
    case = sbt_params.get("case", "-")
    if case == "utube":
        if sbt_params.get("Diameter1") not in [None, "-"]:
            editable_params.append(["Wellbore Diameter Vertical", format_val("Diameter1", sbt_params.get("Diameter1")), "m"])
        if sbt_params.get("Diameter2") not in [None, "-"]:
            editable_params.append(["Wellbore Diameter Lateral", format_val("Diameter2", sbt_params.get("Diameter2")), "m"])
        if sbt_params.get("n_laterals") not in [None, "-"]:
            editable_params.append(["Number of Laterals", format_val("n_laterals", sbt_params.get("n_laterals")), ""])
        if sbt_params.get("lateral_multiplier") not in [None, "-"]:
            editable_params.append(["Lateral Flow Multiplier", format_val("lateral_multiplier", sbt_params.get("lateral_multiplier")), ""])
    elif case == "coaxial":
        if sbt_params.get("Diameter1") not in [None, "-"]:
            editable_params.append(["Wellbore Diameter", format_val("Diameter1", sbt_params.get("Diameter1")), "m"])
        if sbt_params.get("Diameter2") not in [None, "-"]:
            # Diameter2 for coaxial is center pipe diameter
            diameter_val = format_val("Diameter2", sbt_params.get("Diameter2"))
            editable_params.append(["Center Pipe Diameter", diameter_val, "m"])
        if sbt_params.get("center_pipe_thickness") not in [None, "-"]:
            editable_params.append(["Center Pipe Thickness", format_val("center_pipe_thickness", sbt_params.get("center_pipe_thickness")), "m"])
        if sbt_params.get("insulation_thermal_conductivity") not in [None, "-"]:
            editable_params.append(["Insulation Thermal Conductivity", format_val("insulation_thermal_conductivity", sbt_params.get("insulation_thermal_conductivity")), "W/m-°C"])
    
    if sbt_params.get("mesh") not in [None, "-"]:
        mesh_val = format_val("mesh", sbt_params.get("mesh"))
        mesh_labels = {0: "Less Mesh", 1: "Finer Mesh", 2: "Finest Mesh"}
        editable_params.append(["Mesh Fineness", mesh_labels.get(mesh_val, mesh_val), ""])
    if sbt_params.get("accuracy") not in [None, "-"]:
        accuracy_val = format_val("accuracy", sbt_params.get("accuracy"))
        accuracy_labels = {1: "Lower", 2: "Medium", 3: "Higher", 4: "Very High", 5: "Highest"}
        editable_params.append(["Accuracy", accuracy_labels.get(accuracy_val, accuracy_val), ""])
    if sbt_params.get("pipe_roughness") not in [None, "-"]:
        editable_params.append(["Pipe Roughness", format_val("pipe_roughness", sbt_params.get("pipe_roughness")), "µm"])
    
    assumptions.append(["Heat Exchanger Design", case if case != "-" else "-", ""])
    assumptions.append(["Working Fluid", sbt_params.get("fluid", "-"), ""])
    
    if model == "SBT V2.0":
        if sbt_params.get("HyperParam1") not in [None, "-"]:
            assumptions.append(["Inlet Pressure", format_val("HyperParam1", sbt_params.get("HyperParam1")), "MPa"])
    
    if model == "SBT V1.0":
        if sbt_params.get("HyperParam1") not in [None, "-"]:
            hyper1_val = format_val("HyperParam1", sbt_params.get("HyperParam1"))
            hyper1_labels = {0: "Constant", 1: "Variable"}
            assumptions.append(["Mass Flow Rate Mode", hyper1_labels.get(hyper1_val, hyper1_val), ""])
        if sbt_params.get("HyperParam3") not in [None, "-"]:
            hyper3_val = format_val("HyperParam3", sbt_params.get("HyperParam3"))
            hyper3_labels = {0: "Constant", 1: "Variable"}
            assumptions.append(["Injection Temperature Mode", hyper3_labels.get(hyper3_val, hyper3_val), ""])
    
    if sbt_params.get("HyperParam5") not in [None, "-"]:
        assumptions.append(["HyperParameter 5", format_val("HyperParam5", sbt_params.get("HyperParam5")), ""])
    
    if not editable_params and not assumptions:
        return
    
    editable_df = pd.DataFrame(editable_params, columns=["Editable Parameter", "Value", "Unit"]) if editable_params else pd.DataFrame(columns=["Editable Parameter", "Value", "Unit"])
    assumptions_df = pd.DataFrame(assumptions, columns=["Assumption", "Value", "Unit"]) if assumptions else pd.DataFrame(columns=["Assumption", "Value", "Unit"])
    
    if editable_params:
        editable_df.to_excel(writer, sheet_name='SBT Parameters', index=False, startrow=0)
        start_row = len(editable_params) + 3
    else:
        start_row = 0
    
    if assumptions:
        assumptions_df.to_excel(writer, sheet_name='SBT Parameters', index=False, startrow=start_row)
    
    worksheet = writer.sheets['SBT Parameters']
    worksheet.set_column(0, 0, 40)
    worksheet.set_column(1, 1, 20)
    worksheet.set_column(2, 2, 15)

