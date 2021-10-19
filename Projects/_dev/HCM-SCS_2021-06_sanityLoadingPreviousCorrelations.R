
compTable_TTN=
    openxlsx::read.xlsx('/Users/m.wehrens/Documents/Writing/Manuscripts/HCM_SCS/2021_04_Cell_Reports/submitted_manuscript/Dataset_2.xlsx')

compTable_NPPA=
    openxlsx::read.xlsx('/Users/m.wehrens/Documents/Writing/Manuscripts/HCM_SCS/2021_04_Cell_Reports/submitted_manuscript/Dataset_3.xlsx')

toString(strip__chrXX(compTable_NPPA$genes[1:20]))
