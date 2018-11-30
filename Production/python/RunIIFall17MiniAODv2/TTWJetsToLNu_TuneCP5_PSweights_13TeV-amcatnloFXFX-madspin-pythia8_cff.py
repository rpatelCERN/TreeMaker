import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/068FFE22-6743-E811-81AD-0025905B855C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/06A94C8E-6643-E811-9FD6-0025905AA9F0.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0886C1A5-5A43-E811-903A-0CC47A7C3430.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0A20731E-7C42-E811-87A3-0025905A610A.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/1A3F3058-8842-E811-8394-0025905A60AA.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/2E4D2314-9442-E811-9F9F-0025905A60F2.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/327DE03B-5743-E811-B9F9-0CC47A4D76B6.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/32872F96-5343-E811-A628-0CC47A7C3422.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/32D266DF-9742-E811-9E68-0CC47A4D76AC.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/3853DF54-5943-E811-96D9-0CC47A4D763C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/3A86C9A3-6843-E811-9421-0025905A48C0.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/3C0366A5-5A43-E811-A921-0CC47A7C34EE.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/3C9B20BE-5343-E811-9603-0CC47A7C345C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/420F1535-5543-E811-953A-0CC47A4C8EEA.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/46461DFA-5B43-E811-9FC3-0025905B85BA.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/4683D07D-6A43-E811-A85D-0CC47A4D7692.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/4820817C-6743-E811-B144-0CC47A78A3F4.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/487563B4-6143-E811-97C7-0025905B859E.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/5273BF8E-8342-E811-A4AE-0025905AA9F0.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/54F42263-7743-E811-815D-0025905B8580.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/6882B94E-8F43-E811-B975-0CC47A74525A.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/6EFADFAF-5C43-E811-9ABD-0CC47A4D763C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7264C730-7D42-E811-96DE-0025905A6090.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/72FDF71C-6743-E811-95EB-0CC47A7C34C8.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/76987D8C-5543-E811-B01D-0CC47A7C35C8.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7823EC9D-1044-E811-8EB3-0CC47A4C8F06.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/7C764A30-6543-E811-9D2A-0CC47A4D763C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/825330AB-4943-E811-9FFA-0CC47A7C3612.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/82A65A98-5743-E811-811A-0CC47A78A446.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8639D9C6-6343-E811-8187-0CC47A4D7616.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8A4C4BCA-6343-E811-ACD3-0025905A48D8.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8C8B6E10-5A43-E811-92EC-0CC47A4D7692.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8E9C254D-6743-E811-AFB1-0CC47A7C34A6.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/8EC9A7A5-5A43-E811-B3D0-0CC47A4D7628.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/902BF14F-6743-E811-B7EB-0025905B85EC.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/925F9F68-5643-E811-AFA4-0CC47A78A3B4.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/92E45A49-8B42-E811-9AF7-0025905B85AE.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/949C2C87-6643-E811-9FBE-0CC47A7C34EE.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/9C496426-6743-E811-BF4C-0025905B85B8.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/9E27AE1B-4D43-E811-ADC9-003048FFCC16.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/A283C910-5A43-E811-BC8E-0CC47A4D7638.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/A2B8E29A-8C42-E811-A676-0CC47A4D75EE.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/A81D5EAA-8242-E811-BB7D-0CC47A4C8F1C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/AA4E3AAE-5A43-E811-A441-0025905A608E.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/AC787BE7-4A43-E811-9BB1-0CC47A7C34B0.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/ACC29400-6543-E811-BF34-0CC47A7C356A.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B01FF588-6643-E811-9673-0CC47A4D763C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B292219F-6843-E811-AC09-0CC47A4D7604.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B6C31308-6543-E811-8546-0CC47A4C8F2C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/BAABF30F-5A43-E811-BB67-0CC47A4D7670.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/C41E4EB1-9A43-E811-9D97-0025905A608E.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/C8967C0C-5E43-E811-9B6C-0CC47A7C3404.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/D2B388A5-7E42-E811-B141-0CC47A7C3604.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/D849C687-6243-E811-B49D-0CC47A7452D8.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/D867B37A-7343-E811-84A5-0025905A608E.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/DCB7F4A9-8242-E811-A83F-0CC47A4D767A.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/DCC06231-6543-E811-8906-0CC47A4D768C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/EC3B73BE-9643-E811-9863-0CC47A74524E.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/ECBEA291-6643-E811-B451-0025905B8566.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/EEACEA9F-8042-E811-B43E-0CC47A4C8EEA.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/F26E76AF-6943-E811-808A-0CC47A4C8E56.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/FA144220-B242-E811-865C-0CC47A4D76D0.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/FAB198CA-6343-E811-BD20-0025905B85AE.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/FEE713C2-4843-E811-A259-0CC47A4D7638.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/0ABADCD8-C643-E811-98CA-0CC47A7C35D2.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/0AE3E61B-3744-E811-B753-0025905B8586.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/2AC0AD4F-4444-E811-B73B-0025905B8586.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/2E6437B1-B742-E811-9716-0CC47A4D7650.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/329AEEEF-D643-E811-B452-0025905B8582.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/36B849AE-AF43-E811-B328-0CC47A4C8E5E.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/40A1FF61-C143-E811-8824-0CC47A4C8F30.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/48FA897B-C643-E811-8FF0-0CC47A4D7640.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/504F6984-B843-E811-9154-0CC47A4C8E28.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/527F59CD-DE43-E811-BAC3-002618FDA21D.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5E0078E1-D343-E811-A122-0025905B860E.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6449E5AB-BC42-E811-83BC-0025905A610C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6864EA3B-BC43-E811-A138-0025905B8566.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6C81BDDF-C943-E811-BBCB-0CC47A7C35D8.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/7248C462-CD43-E811-A3BE-0CC47A7C34EE.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/76625060-F942-E811-90B9-0CC47A4C8EA8.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/90A9FE81-CC43-E811-9505-0CC47A7C35B2.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/961D3FA6-CF43-E811-A8D0-0CC47A7C3420.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/9EC3DC83-C943-E811-8C1B-0CC47A4D7692.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B01F8B3B-C843-E811-8536-0CC47A4C8F30.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B02C264B-D943-E811-98C6-0CC47A4D760C.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B8BB7501-BF43-E811-A747-0025905B85BA.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C679A89B-DA43-E811-97F4-0025905A6056.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E25A620A-C443-E811-BE86-0025905B857E.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E2679825-B242-E811-8825-0025905A60C6.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F68E5294-C543-E811-A5F4-0025905A60CE.root',
       '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FE7CB2CF-CA43-E811-8050-0CC47A7C34EE.root',
] )
