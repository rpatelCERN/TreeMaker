import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/04A2DC47-1144-E811-88F3-0CC47A4D7600.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/06124D95-D243-E811-A881-0CC47A7C34B0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/1037997B-FE42-E811-B867-0CC47A4D75EC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/141F2ED9-1444-E811-A64B-0025905A60E0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/14CBB9BA-0F44-E811-BDB6-0CC47A4C8F1C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/32A1EB0F-0D44-E811-B3B9-0CC47A4C8E34.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/3EBF2BD6-0543-E811-9211-0025905B85DE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/4AF50071-1B44-E811-950D-0025905A605E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/6015F448-4C44-E811-BFB7-0CC47A4C8F06.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/64E3851C-2444-E811-BF58-0CC47A7C34B0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/74EA11EE-1344-E811-834D-0CC47A4D76CC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/78A37041-0944-E811-A685-0CC47A4C8E64.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/7A960074-0B44-E811-82FA-0CC47A4D75F0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/7E128645-1144-E811-850D-0CC47A4D7630.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/7E925B50-1F44-E811-B74C-0025905B8562.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/8873FDDF-1144-E811-9C90-0025905A6136.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/A80FFBDD-3144-E811-B443-0025905A612A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/BCF183A3-0644-E811-88DA-0CC47A4C8F1C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/BE379460-0A44-E811-AC0B-0CC47A7C3628.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/CCB2368B-0444-E811-9AFD-0CC47A7C35D8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/E0F06F0F-0E44-E811-89D4-0025905B85D2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/F0DA9DAA-1844-E811-BFFE-0CC47A4D7638.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/F6E092FC-0243-E811-A4D9-0CC47A4C8E8A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0441DCFC-DB43-E811-943C-0CC47A4C8E5E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/1A113C99-D643-E811-9D3C-0CC47A7C35F4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/1AA490DD-DC43-E811-8497-0025905A60A6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/1EBC794A-D443-E811-B161-0CC47A78A468.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/32145075-D543-E811-95A8-0025905B85AE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/3893DC01-CE43-E811-BB73-0CC47A7C35A4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/4272126D-CB43-E811-8E12-0CC47A7C3444.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/4602C567-D043-E811-8A4B-0025905A48D8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/561807FB-A743-E811-8F0F-0CC47A7C35D8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/5E424D68-CC43-E811-9D9D-0025905B85CA.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/5E954255-D343-E811-8EE8-0025905B85FE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7203B3E7-F343-E811-98DB-0CC47A4C8EC6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7691044D-3343-E811-92BF-0025905B857E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7885E070-DD43-E811-B1C5-0025905B85FE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/78BB7F59-D443-E811-8B20-0CC47A4D7662.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7C50ADD0-CD43-E811-A903-0CC47A4D7632.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8051284A-D943-E811-9552-0CC47A4C8F30.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/884E68CB-CC43-E811-8DD8-0CC47A4D7616.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8E2A6BC7-D243-E811-99ED-0CC47A4D7662.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/90FA5900-C543-E811-B1A6-0025905B85CC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/9E5B05EB-CE43-E811-BE61-0CC47A4D7692.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/A4AC7F0D-E043-E811-9442-0CC47A4D7640.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/AC2FCDE2-D743-E811-B5F1-0025905B85AE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B49CD603-CA43-E811-B2E2-0CC47A78A41C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/C471DDAC-1644-E811-8CE4-0025905AA9CC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/CE590E06-D243-E811-AD24-0CC47A7C34C4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/D61AFD01-CE43-E811-B1D4-0CC47A7C35A4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/D89CFEDC-CC43-E811-BF3F-0CC47A745298.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/E071550D-E143-E811-B467-002618FDA21D.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/EC7F938C-C843-E811-859B-0CC47A4D7658.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/FA30F016-B943-E811-982F-0025905B8572.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/FE6707EC-E643-E811-A9E8-0CC47A4D766C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/0024ABEF-4E44-E811-BBF9-0025905A60B0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/00FFBF66-8344-E811-8B17-0CC47A4D764A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/02914B21-6D44-E811-8C49-0CC47A78A3EE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/044D4766-4644-E811-B763-0CC47A78A2F6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/164A8A38-5944-E811-800C-0CC47A4D75F0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1834D7B8-6844-E811-B56B-0025905A4964.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1A305AAF-4644-E811-94CD-0CC47A7452D0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1AAA8500-3244-E811-8B31-0025905A612A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/2057B3E9-4E44-E811-AD00-0CC47A78A3B4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/2693F15B-5244-E811-A3FD-0025905B85CC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/26A089E7-4444-E811-A13E-0025905A6134.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/28FBC4D6-5844-E811-B194-0CC47A7C34C4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/36B39F2F-5544-E811-87F9-0CC47A4C8E3C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/38012A2C-4C44-E811-BDFF-0CC47A4C8F26.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/3C4F16E6-4444-E811-B2F6-0025905B85A2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/40972843-4443-E811-B970-0CC47A4C8E5E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/42C739BF-5644-E811-8C2D-0CC47A4C8E34.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/4C3BFE7D-5044-E811-96DE-0CC47A7C3612.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/4C439ECD-5344-E811-A7E3-003048FFCBB2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/4E22FCC5-5344-E811-A563-0CC47A4D7632.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/4ED4A752-5A44-E811-93F1-003048FF9ABC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/52B46CC4-5344-E811-95BC-0CC47A4C8E3C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5876ED59-5244-E811-8EC4-0025905A608C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5A421E8B-5E44-E811-ADFC-0025905A611C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5EA30ECD-5344-E811-9276-0025905A48D8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5EB08D2A-5644-E811-B5AE-0CC47A4D76C6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/62A1740B-3A44-E811-BBDC-0025905B85EC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/666D4459-3D44-E811-814D-0CC47A4D7638.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6A04E723-4944-E811-8E6E-0025905A60F8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/706F14D4-9044-E811-8C88-0CC47A4C8E34.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/70BDB59E-5444-E811-9BFA-0025905A60BC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/7284310D-6044-E811-9BE7-0CC47A4D7698.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/7483429D-4544-E811-8CDF-0025905A60F8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/7C84D7C7-6344-E811-908F-0CC47A7C35B2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/826EBF0F-4E44-E811-939D-0CC47A78A3B4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/8682DF0B-5144-E811-8902-0CC47A7C3424.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/88853C71-4344-E811-AEA7-0CC47A7C35F4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/90DC9E64-4D44-E811-B7EF-0025905B85BC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/9E24F6F2-5744-E811-8BD9-0025905A48D8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A2FBF9D1-5844-E811-9969-0CC47A4D76C6.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A4A69009-5144-E811-8318-0CC47A78A360.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A67686C6-5344-E811-91BA-0CC47A4C8E34.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A86A0265-5C44-E811-8D30-0CC47A7C3444.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A8E1CBC8-4944-E811-BA76-0CC47A7C3410.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/AC43AB0B-4444-E811-B861-0025905A6082.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/ACEFF470-5C44-E811-8328-0025905B858E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B400A279-4244-E811-AE5D-0CC47A4C8F30.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B4D949F4-5A44-E811-BB83-0CC47A4C8E3C.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B810F96D-4C44-E811-9F31-0CC47A7C3612.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B88AC02F-5344-E811-B0C5-0CC47A4D7650.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B8C28FC8-4944-E811-8752-0CC47A4C8ECA.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/BC369C04-1644-E811-AB2E-0CC47A78A468.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/BEFE35A9-5D44-E811-A044-0CC47A4C8F0A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C037543D-4844-E811-9B63-0025905B85FE.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C0D345B0-4744-E811-A747-0025905B8568.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C6916ADC-E644-E811-A595-0CC47A7C35D2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/CE446BCA-4744-E811-95DF-0CC47A7452D0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/CECC0661-6544-E811-8DFE-0025905A60B2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D222F543-3F44-E811-88BF-0025905A6118.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D255E180-4244-E811-B9B0-0025905B85EC.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D45AA3CF-5744-E811-BDA0-003048FFCBB2.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D85C957F-5044-E811-A6F2-0CC47A78A3D8.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E65C1C1A-4944-E811-8267-0CC47A7452D0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E6BDE0A6-4144-E811-830D-0025905A60A0.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E6DC5F1F-6744-E811-8E28-0CC47A4C8F30.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E86C1228-5644-E811-A5A5-0CC47A7C34C4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/EA981C79-4F44-E811-9EA6-0CC47A7C3612.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/EC95220B-3B44-E811-97E1-003048FFD798.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F082EF5F-4B44-E811-B4FF-0CC47A7C3604.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F0B91772-4344-E811-A9D4-0CC47A4D762E.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F0D60A0D-5244-E811-B2C9-0025905B859A.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F463490B-5244-E811-AB76-0025905A60E4.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F4EFEC7B-B644-E811-9DE5-0025905B8592.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F8DB1A66-5944-E811-BE93-0CC47A4C8E34.root',
       '/store/mc/RunIIFall17MiniAODv2/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FE58DC80-4F44-E811-8457-0025905A60F8.root',
] )
