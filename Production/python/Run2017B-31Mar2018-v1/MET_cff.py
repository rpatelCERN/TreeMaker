import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/12CCE365-FD36-E811-B205-008CFA197E8C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/16963797-0937-E811-ABE2-008CFAE45134.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/1AB4F5B1-3637-E811-B23F-008CFAE45300.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/2085229A-1337-E811-830A-008CFA1979AC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/2226ADB0-2E37-E811-90E2-008CFAC93D1C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/261D6CD7-1637-E811-B21A-008CFAE45128.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/36CEFCC1-2A37-E811-8EDE-008CFAC93D94.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/3A2568BB-0637-E811-85ED-008CFAC94300.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/3C0034F7-2B37-E811-AB41-008CFAE45168.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/3E83B61C-FA36-E811-B81A-801844DEDDC8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/566F545C-FF36-E811-8423-008CFAE45040.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/5AEE09DD-2237-E811-8745-008CFAC919EC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/60E6283B-0237-E811-B64D-008CFAE45198.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/686303EF-3437-E811-B399-008CFAE45328.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/6C9390A6-2537-E811-817D-008CFAC94094.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/803643E4-F136-E811-946B-008CFAE453A8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/8265A328-1237-E811-B68F-008CFAC91AC4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/8A2495E7-2B37-E811-80C9-008CFAC91180.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/A29804B6-EE36-E811-B5A4-549F3525B154.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/AC7F5A99-3237-E811-8D97-008CFAE45328.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/ACF73FB4-FB36-E811-8E01-008CFAC93C18.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/B43EB6D0-F636-E811-BB92-008CFAC93C64.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/B4E409C0-1737-E811-A52E-008CFA110C94.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/B879CAF7-F936-E811-AB15-008CFAEBDEEC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/BA7590D2-F136-E811-A080-008CFA197BD0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/C428ACD3-2637-E811-B19E-008CFAE4501C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/CC244FF5-F936-E811-A060-008CFAC91E10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/D010B584-0537-E811-8A60-008CFAC93E14.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/D01DC3A2-0D37-E811-9FA7-008CFAC93DC0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/D404ABDC-F636-E811-89A2-1418774124DE.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/E05A72FA-2737-E811-9996-008CFAC93D24.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/E23B6E40-F436-E811-B07C-008CFAE450B4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/FC4B3573-F836-E811-AE3A-008CFAC915FC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/00273604-AA37-E811-8654-0025905A6080.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/00496277-5737-E811-90B5-008CFA165F58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/00D771BE-AE37-E811-8BB8-008CFAC9429C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/046428B1-7937-E811-81D8-008CFAC94298.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/04A37B49-9A37-E811-B35B-008CFAC93C00.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/0AF85763-DA37-E811-AEA6-008CFAC94044.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/0C46145D-CA37-E811-9E3F-008CFAC91720.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/0CBD980B-6C37-E811-9514-008CFAE452D4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/101B7E58-5537-E811-AC30-008CFAC91CD4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/120E1A7E-B237-E811-BBC1-008CFAC91B58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/1266EE0A-9F37-E811-9E2A-008CFAC91EB4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/12F09421-A337-E811-B0D1-008CFAC9422C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/14BA026E-F437-E811-9905-0025905B861C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/186931ED-B637-E811-AF5F-008CFA11136C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/1A200DE5-B637-E811-B1C1-008CFAC93D64.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/1A57676F-A537-E811-9042-008CFA110C68.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/1ACE9B97-9537-E811-9577-008CFAE45468.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/1C5E6815-7337-E811-848F-008CFAC93FFC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/1C8AD3B2-5737-E811-ACD9-008CFAC93B84.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/1C953225-6C37-E811-B239-008CFA197C10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/1E569F5A-CF37-E811-B2A1-008CFAC93DC0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/22307B5A-CF37-E811-B311-008CFAC93BC8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/22F517EA-BA37-E811-A6EB-008CFA197E0C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/243A921A-9C37-E811-BE1D-008CFAE45440.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/24599584-B037-E811-AFFE-008CFA197E0C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/245D7805-A837-E811-BF39-008CFAEBDBDC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/24A4D1DF-D237-E811-A6DF-008CFAE44CD0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/24C95870-A537-E811-AF32-008CFA110C68.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/28F61FD0-7037-E811-86D1-008CFA5D275C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/2AED06B7-CC37-E811-834B-008CFAE45314.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/2C1CEE99-9537-E811-B47D-008CFAC9402C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/2CF1884F-AC37-E811-9DDE-008CFAC9429C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/2E0F9CAD-8E37-E811-A62D-008CFAE452FC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/300446BE-A037-E811-8108-008CFAC91538.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/321B9722-9637-E811-8A35-008CFAC916B4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/32F027FA-7B37-E811-B7C8-008CFAC91BF0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/380F848B-C837-E811-A3C2-008CFAC94274.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/3E1C377D-9337-E811-9C18-008CFAED6D14.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/3E1E6D65-C337-E811-997F-008CFAE450DC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/3E29146D-B437-E811-B7FF-008CFAE453EC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/4076B008-BF37-E811-A96D-008CFA165F58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/40E31B15-A337-E811-85FE-008CFAE4545C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/428BBE4F-8A37-E811-85DA-008CFAE45440.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/44A4AE51-6E37-E811-A21D-008CFAC94094.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/46078AD6-C637-E811-9FB7-008CFAE45400.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/46374650-AC37-E811-AFDA-008CFAC942AC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/4693D8C5-8E37-E811-8F51-008CFAE45034.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/480A9FB0-A037-E811-B8CF-008CFAC93E24.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/48766F82-6037-E811-995D-008CFA582C28.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/48C6B40F-5C37-E811-B77E-008CFA1979B0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/4A5FABC3-DE37-E811-BAA7-008CFAC91974.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/4ED40F49-9A37-E811-B66D-008CFAC91378.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/50500244-9337-E811-B318-008CFAC910D0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/526EAC07-D337-E811-AE6E-008CFAE45318.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/52A194C5-A937-E811-9A50-008CFAC93D9C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/54BCED6F-CA37-E811-B128-008CFA165F58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/585A1CF0-7E37-E811-8C15-008CFAC93C50.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/58B192AD-8037-E811-A444-008CFAE45314.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/5A54FCBC-D037-E811-B144-008CFA110C94.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/5C030F33-8837-E811-8B58-008CFA197DB8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/5E6A5DAE-A037-E811-BDA9-008CFAE45444.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/6007B1C2-D037-E811-8FE1-008CFAC93EA8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/606FF8B7-CC37-E811-9D58-008CFAC94120.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/6277A404-9837-E811-A1F7-008CFAC919FC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/64628802-BF37-E811-A8FC-008CFAE45024.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/64E4FD9C-6E37-E811-886A-008CFAE453EC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/6645F660-5E37-E811-B6F5-008CFAC93FDC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/66542D8D-C837-E811-8C99-008CFAE45084.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/6E89C568-AE37-E811-AC46-008CFAE45444.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/7037A3EA-BC37-E811-A34E-008CFAC910D0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/705BDCB5-D037-E811-A768-008CFAC919F0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/7240F50B-6C37-E811-80DB-008CFAC91C90.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/743E9FCC-BA37-E811-8D2F-008CFAC93E88.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/78BB76E5-9E37-E811-AAAD-008CFAC940DC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/7A8439C0-CC37-E811-885A-008CFAEBDD90.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/7C53C267-C137-E811-A197-008CFA197C10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/7CEA1B35-D137-E811-A12D-008CFAC941F4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/7E2E18A7-5937-E811-B053-008CFAC93F38.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/7E51C67B-7737-E811-91F3-008CFAC93ED0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/80D2F4E2-B837-E811-B637-008CFAC9429C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/8290F316-A337-E811-80F8-008CFAC91720.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/82DE81CE-8C37-E811-9637-008CFAED6D6C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/82F85A59-AA37-E811-AE0E-008CFA110C68.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/846249AE-A037-E811-9165-008CFAE453DC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/862C751A-8837-E811-922A-008CFA11136C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/8857A19C-CF37-E811-8E42-008CFAC9404C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/88DD6C1A-9C37-E811-AC24-008CFAE452F0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/88E816B1-A037-E811-86D4-008CFAC91204.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/8AF75B58-A837-E811-9DCF-008CFA11136C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/8E34BB02-6337-E811-A496-008CFAE45290.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/9066E3D6-C337-E811-ACF5-008CFAE453EC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/922E6B14-A337-E811-ACE6-008CFAC940D0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/966C4AA3-CF37-E811-9329-008CFAE45134.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/96E7DDB5-CC37-E811-9B01-008CFAC93E88.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/986F47B2-6937-E811-A2FF-008CFAC942AC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/9C811AC1-6237-E811-ACBE-008CFAC93E2C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/9E314681-B037-E811-ACE8-008CFAC942D8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/9E7468AF-CC37-E811-BCB9-008CFAEBDA58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/9E96E914-8D37-E811-A627-008CFAC941A4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/A0BE9C9E-C137-E811-A4FB-008CFAED6D6C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/A403AFC0-F337-E811-AE1B-008CFAC93F74.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/A62F5D76-D537-E811-AF5D-008CFAC94044.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/A64E2E72-D537-E811-9021-008CFAE451BC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/A6609157-DA37-E811-8340-008CFAED6FE8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/AE13F222-A837-E811-9C58-008CFAC941A0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/AECA3185-AE37-E811-9CB2-008CFAE45464.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/B01C3DF9-BC37-E811-944F-008CFAC91490.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/B07CF614-A137-E811-BCA1-008CFAC940DC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/B0A2B649-9A37-E811-A4D3-008CFAC94158.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/B65C1722-6537-E811-8BFC-008CFAE45290.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/B6A7DB02-9837-E811-9DD1-008CFAE453EC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/B823E71A-9C37-E811-9A89-008CFAC93BC8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/BCFCA6DA-9E37-E811-A0F9-008CFAE451BC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/BE4D955E-CF37-E811-AD48-008CFAC91B58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/BE5417A7-8037-E811-BDB9-008CFAC919F8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/C00044CF-8C37-E811-9AC0-008CFAC940BC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/C0CCB134-E337-E811-9AB5-008CFAC94120.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/C48F47A8-5937-E811-9DCF-008CFAEBDB4C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/C4E1F560-C837-E811-B849-008CFAE45124.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/C67A2EEC-9137-E811-BBA2-008CFAC8C264.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/C6DA5549-C637-E811-8F91-008CFAED6D6C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/C8156D7B-7737-E811-878E-008CFAC94234.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/C8445B20-9C37-E811-9020-0025905A48D0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/CA11BE6C-C337-E811-815F-008CFAC93CB4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/CCB7DFED-C337-E811-AF39-008CFAE450F0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/CCC88AE7-B837-E811-A1F4-008CFA11136C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/CEA746AE-A037-E811-8074-008CFAC91DD0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D0128B5A-CF37-E811-968E-008CFAC91D8C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D2440346-E737-E811-921D-008CFAC9414C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D2528F61-A537-E811-818C-008CFAC93E18.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D2671A16-8837-E811-8D91-008CFA1979EC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D432DEB3-6837-E811-90EC-008CFA197C10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D43E16DA-9E37-E811-BAC4-008CFAE453A4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D46F8575-C237-E811-B93B-008CFAC91B58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D61F2306-D337-E811-948D-008CFAE452A0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D647F4E6-9E37-E811-BB34-008CFAC942AC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D8AB9406-D337-E811-B5AF-008CFAC93B38.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D8C44C32-E537-E811-BFAE-008CFAC93BB0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/D8D9110E-A837-E811-B3AC-008CFA197B74.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/DCCDF4FE-9037-E811-BECC-008CFAC93FFC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/DED31276-6737-E811-BF22-008CFA197C10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/E0635272-D537-E811-ACA0-008CFAC93CD8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/E0A2F658-AC37-E811-9A7B-008CFAE4507C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/E0F968EB-C337-E811-8882-008CFAED6D6C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/E2B44284-C137-E811-A4A8-008CFAE45034.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/E635C3C0-A737-E811-B403-008CFAEBDB6C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/E6FF5282-B437-E811-8A47-0025905A60C6.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/EC05FDCF-8337-E811-9329-008CFAE44CD0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/EC262681-B037-E811-A92D-008CFAE45058.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/EE615EBD-A937-E811-B559-008CFAC919F0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/EEBC2F66-CF37-E811-BE46-008CFAC910D0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/F246CCF8-BC37-E811-943E-008CFAC93FD8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/F26E3265-AE37-E811-B9A4-008CFAE45058.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/F2EB617D-B237-E811-B110-008CFAC918F0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/F4F49516-B737-E811-A472-008CFAC941A0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/F8EB3164-A537-E811-A7AA-008CFAC94030.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/FCD9F09C-CA37-E811-9E73-008CFAEBDA58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/80000/FE541415-A337-E811-83A9-008CFAC93E18.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/0084732F-3F37-E811-8A86-008CFAC8DB40.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/02100337-1F37-E811-B701-008CFAE45278.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/02957166-7437-E811-9776-008CFAC9411C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/04AACA1F-2837-E811-9498-008CFAC93EC0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/04B41A11-3737-E811-AA79-008CFAC91BF0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/06D4B81A-4A37-E811-AC9F-008CFAC93C54.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/080D994F-3137-E811-BE7C-008CFAC93BD0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/0AC03A9B-1137-E811-8F14-008CFAC93BC4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/0C422D7F-0237-E811-BCBC-008CFA197C10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/0C430F74-0237-E811-BB28-008CFAC93C64.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/10CE1BB4-1737-E811-9253-008CFAC940BC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/10F55C23-2C37-E811-8353-008CFAC91BF0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/12119F1B-3D37-E811-B707-008CFAC93B94.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/1216BF3B-1B37-E811-B32F-008CFAE4515C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/1257D35F-4C37-E811-9B2C-008CFAC93C54.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/1650F834-4137-E811-BAC1-008CFAC8DB40.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/16724C6F-3137-E811-A630-008CFA197E0C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/186221B3-1737-E811-A72E-008CFAC916E4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/1C6CB412-3037-E811-AD1E-008CFA197BD0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/1CAB90EC-3A37-E811-BBAE-008CFAE45078.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/1CC7358D-1537-E811-8BC6-008CFAC93C9C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/2294E7FE-5A37-E811-B5A4-008CFAC91B94.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/22A95A4A-0A37-E811-8595-008CFAE4535C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/22C0EB76-5337-E811-9547-008CFAE45214.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/242672FF-2C37-E811-8EC2-008CFAE45198.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/247D0C84-F836-E811-8D57-008CFAC93B60.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/2487F7EB-1837-E811-B1BF-008CFAE4507C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/2659F31B-3737-E811-B145-008CFAE451E0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/2AA6326B-6137-E811-BB0D-008CFAC94154.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/2AE021B4-FE36-E811-8D1E-008CFAC93E8C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/2C7DFFA1-5737-E811-B845-008CFAE45348.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/2E828B4B-2A37-E811-BAEE-008CFAC91940.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/2EA433B2-1737-E811-A42F-008CFAC94274.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/32A63CE4-0D37-E811-BB4D-008CFAE45078.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/380E0B77-0237-E811-86C3-008CFAEBDD90.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/38558E54-0A37-E811-9159-008CFA1C9308.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/3E4A2519-3B37-E811-A2EB-008CFAC91A00.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/40F4D234-F836-E811-88D2-008CFAED6D14.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/42884553-3F37-E811-A6F6-008CFAC93B8C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/44264D96-1337-E811-88BF-008CFAE45340.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/446BFCF7-1737-E811-97BB-008CFAC94134.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/468F5639-1F37-E811-A286-008CFAE45280.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/46A91DB0-1737-E811-8163-008CFAE45368.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/4857D693-4337-E811-93AA-008CFAC8DB40.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/4A38966F-4C37-E811-9B9E-008CFA165F58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/4A393697-1537-E811-A582-008CFA197B74.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/4C28A1EE-3437-E811-A95A-008CFAE45328.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/4EB90914-3237-E811-8D18-008CFAE45070.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/4EB9FBE7-0237-E811-8489-008CFAC93F58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/528C2D4E-1B37-E811-AE3D-1866DAEA79C8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/541B059D-1337-E811-9E9B-008CFAE45080.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/56C41203-6337-E811-9E6D-008CFAC94080.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/5A0B202A-1B37-E811-92CE-008CFA1C9308.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/5E9738C5-3637-E811-A46F-008CFAE450D4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/60236456-1937-E811-9426-008CFAC9406C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/6048EF1A-4A37-E811-9CE4-008CFAC941A4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/623FB57D-2837-E811-9147-008CFAC91E10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/649BC51A-4A37-E811-B2FD-008CFAC91DCC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/66A596FF-3C37-E811-B549-008CFAC93DA8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/68F4B833-5E37-E811-BD94-008CFAC94064.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/6A4E90F8-0B37-E811-BC5C-008CFAC93DD4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/6A63D23B-2137-E811-8CE5-008CFAE45280.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/6CCAA6EE-3537-E811-8DDE-008CFAC94080.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/70150C8F-1D37-E811-9834-008CFA1979AC.root',
] )
readFiles.extend( [
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/7097C80A-5437-E811-9F4B-008CFA197C10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/70B4DAAD-3637-E811-8125-008CFAE45290.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/7205804E-2A37-E811-935D-008CFAE452FC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/726BBABF-3637-E811-9BA2-008CFAE451D8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/72D4E79E-1137-E811-9F7D-008CFAC94274.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/74863851-0A37-E811-B715-008CFAE45328.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/76410A34-F836-E811-B356-008CFAC91E10.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/78354CA6-7F37-E811-A8CC-0CC47A4D760C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/78646527-1B37-E811-BCAF-008CFAE4507C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/78C7F18C-9337-E811-B130-0CC47A4D769A.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/78ED5D38-2E37-E811-80A8-008CFAC93DC8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/7C9FF57B-3837-E811-B4AA-008CFAC91ED8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/7EDE85C5-1337-E811-883F-008CFAE4531C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/8088BDD7-3637-E811-A394-008CFAC93CD8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/8437F196-1837-E811-ABE6-008CFAC941C8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/848366DF-4537-E811-8EA4-008CFAC94138.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/84845598-0637-E811-9ADF-008CFAE45094.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/8493C12C-FC36-E811-A55B-008CFAC93E8C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/860A635B-2637-E811-A1A4-008CFAC93EDC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/86BDCC3C-1937-E811-AC0B-008CFA1979AC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/887A4A35-F836-E811-B418-008CFAE45048.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/8899679E-1137-E811-B638-008CFAC93E8C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/8CFDCA59-2637-E811-8F9F-008CFAC940A0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/8E1A678C-1537-E811-ABB8-008CFAC93C9C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/8E2F9E21-3937-E811-A3B9-008CFA197A90.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/8E8ACE25-2C37-E811-9278-008CFAC93FD4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/90F022DF-6037-E811-947C-008CFAC94208.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/928B0345-6937-E811-B533-008CFAC94208.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/92A7BE0B-6C37-E811-807B-008CFAC93E88.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/949AA696-1337-E811-A169-008CFAE45340.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/94B65B46-2E37-E811-9AB9-008CFA197AE4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/961476AA-4A37-E811-9A45-008CFAC94288.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/965D43C7-1737-E811-AC7B-008CFAE453E8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/96663935-1937-E811-86ED-008CFAC8DB40.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/989150E4-3437-E811-B43F-008CFAC91BF0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/9E53BD58-3137-E811-AD85-008CFA197AF4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/9E67B3E7-2437-E811-A03F-008CFAE45200.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A03167BD-1537-E811-8B97-008CFAEA1710.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A07C478D-0637-E811-AA03-008CFAC941A0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A08AEE2D-FC36-E811-800B-008CFAE45070.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A25B85BE-2E37-E811-A664-008CFAC93FD4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A451EA33-F836-E811-8407-008CFAC93DC0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A49F69CE-4E37-E811-9C70-008CFAE45278.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A4C1790E-5A37-E811-9568-008CFAC91B58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A4EF7FD3-5B37-E811-A784-008CFAC916F0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A8047666-0A37-E811-9F8B-008CFAC93D4C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/A82D570E-3037-E811-BB5A-008CFAE4526C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/AA8CDB3B-3B37-E811-A09B-008CFAC916E4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/AC5C3545-FC36-E811-98B8-008CFAC93DD4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/AC875DFB-0B37-E811-A98C-008CFA197980.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/AC94924E-3037-E811-8AD0-008CFA11136C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/B2023937-1937-E811-92CE-008CFAC8DB04.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/B2230AB5-1737-E811-876D-008CFA1C9308.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/B22E6DD1-0A37-E811-980D-008CFAE45078.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/B2BBEF99-5537-E811-B6B1-008CFAC94298.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/B8949847-1B37-E811-8C3B-008CFAE45040.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/B8A10538-1F37-E811-AB6E-008CFAC93F08.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/BCCEFECD-1D37-E811-AD5D-008CFAE45070.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/BCEBCB8B-7F37-E811-9F9A-008CFAC93FFC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/BE756E9E-1137-E811-9DB6-008CFAE45368.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/BE8291D3-5B37-E811-8E5F-008CFAC916F0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/C07337C7-2437-E811-8159-008CFAC93F9C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/C63A74EA-1537-E811-9E6C-008CFA197BDC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/C8FDC99C-1137-E811-B969-008CFAC93E08.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/CA48E233-FC36-E811-81F0-008CFAC93E48.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/CA4E604A-3437-E811-8EF2-008CFAE45070.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/CCBB6B48-1B37-E811-991D-008CFAC915FC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/CE021F56-1B37-E811-908A-782BCB539B14.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/CE7D5589-1D37-E811-9655-008CFAC915FC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/D06A7EA0-1137-E811-9BED-008CFA197980.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/D0CD0263-7637-E811-9F2C-008CFAC915FC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/D0F31909-5C37-E811-A4B0-008CFAC94018.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/D4066287-3837-E811-8F1F-008CFA197DC8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/D45AF31B-3737-E811-BE18-008CFAE451E0.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/D65AA24F-2A37-E811-B915-008CFAE450B4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/D830E3A0-0637-E811-B3F9-008CFA197AF4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/DA684C9D-1137-E811-A8F7-008CFAE4531C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/DA85F18B-1537-E811-8461-008CFAC93D90.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/DC99AF96-1337-E811-897E-008CFAE45340.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/DE50297B-1F37-E811-AC93-008CFAC93B5C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/E43F048B-FE36-E811-9B9C-008CFAC93E48.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/E4A64924-2C37-E811-96C1-008CFAE4502C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/E4DFA98D-1537-E811-9211-008CFAC93C9C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/E6526B54-3137-E811-8BA0-008CFA165FE4.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/E6A1A973-0237-E811-8AB3-008CFAE4535C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/ECF6E74D-2A37-E811-AED8-008CFAC93B5C.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/EECE8B96-1337-E811-BC75-008CFAE45340.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F0B05736-1937-E811-BCBC-008CFAE45134.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F0CCB39F-5537-E811-AC1D-008CFAE45348.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F0E87634-1937-E811-9FF1-008CFAE45410.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F26D7C72-FE36-E811-BD0B-008CFAC91A38.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F46FF0B0-2437-E811-AE56-008CFAC94014.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F4C0EF2B-0E37-E811-A818-008CFAC93F58.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F4C7324A-0A37-E811-BCE6-008CFAE451DC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F4E06643-2837-E811-8DE0-008CFAC91C90.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F808A69C-1337-E811-8B95-008CFAC91ED8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F80E358D-0637-E811-B939-008CFAE450CC.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/F8F3B9F6-5537-E811-9917-008CFAC93E68.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/FA8D52E2-0D37-E811-BAA5-008CFAC93DB8.root',
       '/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/90000/FE13E873-0237-E811-ACE8-008CFAE4528C.root',
] )
