import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/0CE70EF2-7042-E811-AFCE-0242AC1C0503.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/3441323F-7442-E811-A8C1-0242AC1C0501.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/3E6E6D85-7042-E811-8567-0242AC1C0501.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/44C05F64-7942-E811-ABB9-0242AC1C0504.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/4681F65D-7142-E811-A2C4-0242AC1C0503.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/4EC4619F-7442-E811-A874-0242AC1C0500.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/6E55B71F-7542-E811-A2C9-0242AC1C0503.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/7E341501-7E42-E811-9E2F-0242AC1C0504.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/80660F95-6B42-E811-BDA5-0242AC1C0500.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/8C091CB6-7F42-E811-8F88-0242AC1C0501.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/94BD23BC-4242-E811-9D72-0242AC1C0500.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/A623FC90-7242-E811-A1A1-0242AC1C0501.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/AC5EC697-7442-E811-94E6-0242AC1C0500.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/C2A77521-7842-E811-A79C-0242AC1C0500.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/C4760E90-7742-E811-847E-0242AC1C0502.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/DEA11B9D-FA41-E811-ACD0-0242AC1C0502.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/F6192C58-7642-E811-A820-0242AC1C0501.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/F6F88108-7342-E811-A7B5-0242AC1C0503.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/FAAA4E21-7842-E811-BE4B-0242AC1C0500.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/0642FF77-7243-E811-9EC5-44A842CFC98B.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/0699705E-7343-E811-9893-484D7E8DF09F.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/087D997D-7243-E811-B47C-484D7E8DF0E0.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/0EFD1944-7343-E811-8971-484D7E8DF107.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/20217077-7243-E811-A04B-44A842CFC98B.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/300C88F5-7243-E811-9A42-B499BAAB50A0.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/30C463A7-7343-E811-928D-6C3BE5B59218.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/36766F43-7343-E811-AAAC-484D7E8DF0FA.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/3C778035-7343-E811-A9AA-484D7E8DF0B9.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/3E5ACF1F-7343-E811-8EBE-484D7E8DF0B9.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/42C8F979-7243-E811-ADF1-44A842CF05B2.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/46E0A5CF-7243-E811-B4FE-44A842CFCA27.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/484DF54C-7343-E811-8B6D-B499BAAC00E0.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/562D0B81-7243-E811-8BC1-001F290789D6.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/60580B23-7343-E811-B042-44A842CF05E6.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/68D24124-7343-E811-89C7-44A842CF05E6.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/6A6A427B-7243-E811-B18B-484D7E8DF0AC.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/6AA3FCB0-2544-E811-8490-0242AC1C0503.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/6C440B81-7243-E811-B84E-001F290789D6.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/7454DE86-7243-E811-9724-6C3BE5B50178.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/764254FD-7243-E811-8286-484D7E8DF092.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/7697FD4E-7343-E811-B9F2-B499BAAC0A22.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/7A07E882-7243-E811-B9EA-6C3BE5B59058.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/82373D22-7343-E811-AF62-44A842CF058B.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/864804E7-7243-E811-9EB6-44A842CFC98B.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/8ACF6EFB-7243-E811-81CB-001F29087EE8.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/8EEC067A-7243-E811-AF42-44A842CF05B2.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/A86CB8D3-7243-E811-BCD5-484D7E8DF06B.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/AA26777C-7243-E811-BB39-484D7E8DF114.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/AC04E4FB-D243-E811-87DA-0242AC1C0501.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/B2BCFA4A-7343-E811-9884-44A842CF05A5.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/C4F19109-9043-E811-B8E9-484D7E8DF0AC.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/CC29983F-BE43-E811-AE3F-0242AC1C0501.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/DABDDE74-7343-E811-890B-001F2908CE36.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/DE549331-1144-E811-984B-44A842CFD667.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/E0D220C7-7343-E811-892A-B499BAAC0A22.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/E2FA73F0-7243-E811-8FE1-6C3BE5B50178.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/E4A97D86-7243-E811-A105-6C3BE5B50178.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/EEE3DDE1-7243-E811-9DC8-B499BAAC055E.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/F432FA4C-7343-E811-A39F-B499BAAC00E0.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/F620777C-7243-E811-ABEC-484D7E8DF114.root',
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/FA41FCEF-7243-E811-99C0-B499BAAC055E.root',
] )
