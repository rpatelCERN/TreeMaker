import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/1E4C4F31-7374-E511-875C-0025905A6138.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/24BBA38D-F774-E511-9A1A-0025905B85D6.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/34BE5917-7374-E511-B279-0025905A612C.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/4AA0FEBB-7274-E511-928B-0025905B858A.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/50D108BC-7274-E511-A1BB-0025905B858A.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/56008BFB-6774-E511-94E3-0025905A60D2.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/646E27EC-7274-E511-B3CF-002590593902.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/727BCBC8-7274-E511-85A4-00261894386E.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/7C991971-6874-E511-93B3-0025905A6066.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/82C67545-7374-E511-A668-0026189438CC.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/8C394E37-6874-E511-A9E3-0025905B85F6.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/B804389C-7274-E511-8F91-0025905A60EE.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/BE6D2B01-6F74-E511-9F4F-0026189438F7.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/C26C8E46-7374-E511-8381-0025905A608C.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/E2951C47-7374-E511-A42C-002618943956.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/FA0F81D8-7274-E511-9A2F-0025905B858C.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/0EBF04E8-7174-E511-8E21-0025905A60BE.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/123B4A40-7274-E511-A9A0-0025905A48F0.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/14254C00-7474-E511-B79A-0026189438F4.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/185A3501-7474-E511-8C61-002618943947.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/206CCB3E-7274-E511-A977-0025905A6070.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/3A66B13F-7274-E511-9AC7-0025905A6094.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/44BD8244-7274-E511-858D-0026189438F4.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/44C83043-7474-E511-BF93-002618943807.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/48C47F17-7474-E511-A0C3-00261894395B.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/50982B3A-7474-E511-99D3-0026189438F2.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/50FB353D-7274-E511-AB5D-002618943940.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/52E87138-7474-E511-9661-0025905A609E.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/5434563F-7274-E511-9804-0025905B8572.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/58E30E42-7274-E511-A8E8-0026189438D8.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/5C8C51F1-7374-E511-B216-0026189437EB.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/5E23F608-7474-E511-BCC6-0026189438D8.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/6631C1EB-7374-E511-8ADE-0025905A48D0.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/68E38F42-7274-E511-88E9-0025905A60D0.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/6A54A13B-7274-E511-B61E-00261894395B.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/6CC0FC3F-7274-E511-ADE1-0025905A48EC.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/8EEA563F-7274-E511-9A9F-0026189438E4.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/92CFC13E-7474-E511-9246-0025905964BE.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/96570445-7274-E511-84C1-002618943947.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/989B3E24-7274-E511-AC03-0025905A60E4.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/9CDB5E42-7274-E511-AD96-0025905A612A.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/9E4EB53F-7274-E511-A4EE-0025905A60BE.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/AAD3F644-7274-E511-A0C9-003048FFD770.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/B0BE1C36-7474-E511-8584-0026189437EC.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/B27C6041-7474-E511-AAB7-003048FFD736.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/B4E41D3F-7274-E511-BF17-002618943986.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/C0892CD0-7174-E511-9DF1-0025905A60A0.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/C0F9433E-7274-E511-B2B0-0025905A60DE.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/C4EC293E-7274-E511-9A23-0025905A6092.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/DECC123D-7274-E511-B0D1-0026189437EB.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/E2364423-7474-E511-A140-0025905A6094.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/E8CF78EB-7374-E511-A0A9-002618943986.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/EA062DE0-7374-E511-A248-002618943865.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/EC59543D-7274-E511-BCD2-002618943865.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/ECEA5C06-7474-E511-8344-0026189438E4.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/F0E3E941-7274-E511-9C03-0025905B858C.root',
       '/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/80000/FEFD7047-7274-E511-BD6D-003048FFD7C2.root' ] );


secFiles.extend( [
               ] )
