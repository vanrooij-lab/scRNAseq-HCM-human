
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_sessions/HCM_2020may_re-analysis-v4.Rdata')

groupedSCS_SinglePatients = list(
    patient1Mod=groupedSCS$patient1Mod,
    patient2Mod=groupedSCS$patient2Mod,
    patient3Mod=groupedSCS$patient3Mod,
    patient4Mod=groupedSCS$patient4Mod,
    patient5Mod=groupedSCS$patient5Mod)

save(list='groupedSCS_SinglePatients', file='/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Rdata_previous/groupedSCS_SinglePatients.Rdata')

groupedSCS_PatientAllMod = list(patientAllMod=groupedSCS$patientAllMod)
save(list='groupedSCS_PatientAllMod', file='/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Rdata_previous/groupedSCS_PatientAllMod.Rdata')
