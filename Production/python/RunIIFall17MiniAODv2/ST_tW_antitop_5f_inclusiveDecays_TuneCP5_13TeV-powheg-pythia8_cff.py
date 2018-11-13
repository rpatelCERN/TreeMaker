import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/04B402CD-FA42-E811-80FD-AC1F6B1AF148.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0C5C3506-0743-E811-900F-A0369F83630C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/10AAB630-0943-E811-A1B3-A0369F83630C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/18E6318D-0143-E811-999A-A0369F836288.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/2822DAF8-F642-E811-9BD7-1CC1DE1CDD20.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/2ADC2266-BB43-E811-BD1C-00266CFF0AF4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/420BBEB0-0243-E811-A320-7CD30AC03712.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/48D37F50-FB42-E811-8049-008CFAF554D2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/54AE12CE-FE42-E811-8048-6CC2173BBE90.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/5A04C5A1-0343-E811-BF91-C4346BBC9BB0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/6076FF20-FF42-E811-86CD-A0369F83630C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/60DB8807-FD42-E811-B06A-A0369F8363DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/6C7D8CD1-FE42-E811-A514-A0369F836280.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7067E08C-0143-E811-9F96-A0369F836342.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/784A0271-FC42-E811-BD11-A0369F83637E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7C436F1D-FF42-E811-B307-A0369F836342.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7CBED001-FD42-E811-8ADC-A0369F8363DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/806DA5F8-0043-E811-8164-008CFAF554D2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8412B49E-FD42-E811-B38F-C4346BBC1498.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8E0CDA8E-0143-E811-A7C2-AC1F6B1AF03C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/9CE2A440-0543-E811-AF36-C4346BBC9BB0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/A2E15E6E-FA42-E811-A8FE-7CD30ACE2445.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B20BD89E-A342-E811-921A-A0369F83628A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/BA96C863-0043-E811-9D18-AC1F6B1AF148.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/C870518D-0143-E811-B289-A0369F83630C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/CEAE29AF-0943-E811-9C1F-A0369F8363DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/EA5AB3A1-A742-E811-AD1C-AC1F6B1AF1CE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/FE132722-FF42-E811-8F23-008CFAF5550C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/000DEAF4-8B43-E811-9784-1CC1DE1D03DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/168ED6E7-8843-E811-840E-C4346BC7EDD8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/1CBF8020-7E43-E811-9CEC-C4346BBC1498.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/248AC1BC-8243-E811-BB6A-A0369F83642A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/28A3568F-8043-E811-9203-C4346BBC9BB0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2AD32357-8743-E811-8E14-008CFAF5550C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/303E3824-8543-E811-BFDE-C4346BBC1498.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/34696FA2-1E43-E811-94F2-A0369F836364.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4E8B36B2-8D43-E811-A838-7CD30AC03722.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/6670F2AB-8A43-E811-B1F9-7CD30AC03722.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/6A1D2C68-8E43-E811-8C50-1CC1DE056008.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/722002B7-9543-E811-A7C4-A0369F8362E8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7692B19C-9943-E811-84F1-7CD30ACE2445.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7E64A5D9-8943-E811-AD0E-00266CFE8A04.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/88072ABE-8243-E811-ABFA-C4346BBC9BB0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/8C5EBF04-8843-E811-94C7-AC1F6B1AF194.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/98F231F2-8C43-E811-90E7-A0369F8362E8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A030D649-8C43-E811-AE72-7CD30AC03722.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A2791133-8143-E811-B2AD-00266CFFC7E0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/AA396C6A-7C43-E811-A362-C4346BBC1498.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/B26E12B2-8343-E811-B8E6-AC1F6B1AF1CE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/B423DE81-8643-E811-B68B-AC162DACC328.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/B4CE69C2-8243-E811-81C6-AC1F6B1AF148.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/BC516CFD-8743-E811-B8DB-C4346BBC1498.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/CC8B9AD2-8543-E811-AF69-1CC1DE1D03DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D64078DA-9043-E811-9A0B-00266CFFC940.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F459FE45-0443-E811-A31B-A0369F8363DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F4E4C0EF-7F43-E811-8DB9-008CFAEEAD4C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F8BB34A0-FE42-E811-9CE6-A0369F83637E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/021E1DE3-1044-E811-8449-AC1F6B1AF144.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/04EE4B2A-7844-E811-A043-008CFAEEABF8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/2251040A-9B44-E811-B725-7CD30AC036FE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/280ACFD4-D244-E811-A565-008CFAEEABF8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/3298F42B-7844-E811-A47F-7CD30ACE1239.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/4A1AD81B-1544-E811-9766-001E67E5E8B6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/4C0A6415-7544-E811-A429-008CFAEEABF8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5899719D-7644-E811-BC5C-00266CFEFE1C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/900C1DF4-0A44-E811-935C-7CD30AC0370E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/9E7555C0-7544-E811-80BF-008CFAF5550C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FA8665CA-7A44-E811-BD71-AC1F6B1AF18A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FC9457A7-7344-E811-8702-AC1F6B1AF01A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/04D3B302-4543-E811-A5E8-7CD30ACE2445.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/060CF388-4E43-E811-8E31-C4346BC80410.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/0C9C9877-6043-E811-8950-00266CF65AC4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/109086F3-4943-E811-BB56-7CD30ACE1239.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/10AF6D11-5D43-E811-9D60-1CC1DE1D03DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/10F05A4D-6143-E811-94D0-AC1F6B1AF186.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/12AD4A53-5743-E811-99AF-AC1F6B1AEF94.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/1813F4E3-5643-E811-AFD2-1CC1DE048FD0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/18E27537-5D43-E811-9041-00266CFE8A04.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/1A033669-5743-E811-806A-AC1F6B1AEFEE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/1C7FFB4D-9942-E811-AE38-A0369F8362E6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/24DDC09A-4843-E811-AE7F-A0369F836342.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/285C82C9-5C43-E811-B93A-AC162DA87230.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/2A05D0D9-6243-E811-A4C4-00266CFFC7E0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/2C6E36B4-5243-E811-B3A2-A0369F836342.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/2EE7DED6-6143-E811-B6DA-001E67DBE26E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/30976B87-5443-E811-BFF6-00266CFFC7E0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/321E93C5-5643-E811-956B-AC1F6B1AEFFC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/361CA3E2-4543-E811-8B39-1CC1DE056008.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/3E171C6A-5143-E811-95B8-78E7D1E49636.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/40817721-5543-E811-9F72-C4346BBC1498.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/48E598F9-9042-E811-BA39-A0369F836280.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/4A061B2C-5E43-E811-BCCF-A0369F836342.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/4ABF1BD1-6D43-E811-9C8F-AC1F6B1AF148.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/4CEEB684-B742-E811-B8CE-C4346BC7EDD8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/52B023DD-4C43-E811-AAAC-00266CFFBF88.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/580DFC48-4B43-E811-9A34-00266CFFBCD0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/5A27853D-8F42-E811-92CD-C4346BBC1498.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/60E2B945-4043-E811-B793-A0369F836280.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/62786C6A-4443-E811-BFCA-1CC1DE056008.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/64B92922-5043-E811-8B30-1CC1DE1D03DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/6A0FD043-4043-E811-9616-A0369F83642A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/7CCD5662-A342-E811-8E46-AC162DACB208.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/84DE0731-4143-E811-8D47-7CD30ACE2445.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/88E8C060-4A43-E811-9FE5-AC1F6B1AF01A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/8E7C72F2-4443-E811-914D-00266CFFC980.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/9024CB81-5343-E811-94F7-A0369F836280.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/909EA033-4843-E811-B6FB-AC1F6B1AF1CE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/98113ACC-5443-E811-8E5E-00266CFEFF04.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/982BEF64-9342-E811-8D3D-AC1F6B1AF03C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/9A5F227E-3D43-E811-9C59-A0369F836280.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/AE8C4EDE-4243-E811-A88F-C4346BC70EC8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B220657B-3D43-E811-9F9F-AC1F6B1AF148.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B2E2A418-6743-E811-9FEF-A0369F8363DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B4BAAFE7-8B42-E811-84F2-7CD30AC03712.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/BEE3973A-3643-E811-B847-7CD30ACE2445.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/C4F68E24-5043-E811-B32D-A0369F836280.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/C8A0A677-5A43-E811-A327-00266CFFBCD0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/CA8DD508-5643-E811-B8A2-A0369F83639C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D02E3662-4A43-E811-8C45-A0369F8362E6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D2EA6D60-4243-E811-B57A-AC1F6B1AF148.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/DA168C8C-4643-E811-8023-A0369F836342.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/DC2112EF-7343-E811-98BE-00266CFFCCBC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E29D17BB-4343-E811-8721-00266CFEFF04.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E62A99A4-A742-E811-87A3-7CD30AC03712.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/EA1CE572-5A43-E811-B413-001E67E5E8B6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/EAD18029-7643-E811-8AE3-AC1F6B1AF18A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/EC767D6D-4D43-E811-902F-AC1F6B1AF01A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/ECCD432E-9842-E811-A88A-1CC1DE1CDD20.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/EE04B269-A642-E811-9DDA-AC1F6B1AF140.root',
] )
