#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import pandas as pd

def write_excelsheet(df_summary, df_subsurf_res_mass, df_subsurf_res_time, df_econ, geoCLUSTER_results_pathname):

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

    writer.close()

