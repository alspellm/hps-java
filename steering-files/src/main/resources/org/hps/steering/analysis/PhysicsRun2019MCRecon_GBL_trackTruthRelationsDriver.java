<?xml version="1.0" encoding="UTF-8"?>
<lcsim xmlns:xs="http://www.w3.org/2001/XMLSchema-instance" xs:noNamespaceSchemaLocation="http://www.lcsim.org/schemas/lcsim/1.0/lcsim.xsd">
    <!-- 
      @brief Steering file that will be used for the reconstruction of 2019 MC. 
      @author <a href="mailto:omoreno1@ucsc.edu">Omar Moreno</a>
    -->
    <execute>
        <driver name="EventMarkerDriver"/>
        <!-- Ecal reconstruction drivers -->        
        <driver name="EcalRunningPedestal"/>
        <driver name="EcalRawConverter" />
        <driver name="EcalTimeCorrection"/>
        <driver name="ReconClusterer" />
        <driver name="CopyCluster" />
        <driver name="HodoRunningPedestal"/>
        <driver name="HodoRawConverter"/>
        <!-- SVT hit reconstruction drivers -->
        <driver name="RawTrackerHitSensorSetup"/>
        <driver name="RawTrackerHitFitterDriver" />
        <driver name="TrackerHitDriver"/>
        <driver name="HelicalTrackHitDriver"/>
        <!-- Track finding and fitting using seed tracker. --> 
        <driver name="TrackReconSeed123Conf4Extd56"/>    
        <driver name="TrackReconSeed123Conf5Extd46"/>
        <driver name="TrackReconSeed567Conf4Extd123"/> 
        <driver name="TrackReconSeed456Conf3Extd127"/> 
        <driver name="TrackReconSeed356Conf7Extd124"/>  
        <driver name="TrackReconSeed235Conf6Extd147"/>  
        <!-- <driver name="TrackReconSeed234Conf6Extd157"/> -->
        <driver name="MergeTrackCollections"/>
        <driver name="GBLRefitterDriver" />
        <driver name="TrackDataDriver" />
        <driver name="KalmanPatRecDriver"/> 
        <driver name="TrackTruthMatching_KF" />
        <!-- <driver name="TrackTruthMatching_GBL" /> -->
        <driver name="ReconParticleDriver" /> 
        <!-- <driver name="ReconParticleDriver_Kalman" /> -->
        <!-- <driver name="TrackClusterTruthMatchingDriver" /> -->
        <driver name="TrackTruthRelationsDriver"/>
        <driver name="LCIOWriter"/>
        <driver name="CleanupDriver"/>
    </execute>    
    <drivers>
        <driver name="TrackTruthRelationsDriver" type="org.hps.recon.utils.TrackTruthRelationsDriver">
            <enablePlots>True</enablePlots>
        </driver>

        <driver name="EventMarkerDriver" type="org.lcsim.job.EventMarkerDriver">
            <eventInterval>10</eventInterval>
        </driver>
        <!-- Ecal reconstruction drivers -->
        <driver name="EcalRunningPedestal" type="org.hps.recon.ecal.EcalRunningPedestalDriver">
            <logLevel>CONFIG</logLevel>
        </driver>
        <driver name="EcalRawConverter" type="org.hps.recon.ecal.EcalRawConverter2Driver"/>
        <driver name="EcalTimeCorrection" type="org.hps.recon.ecal.EcalTimeCorrectionDriver"/>
        <driver name="ReconClusterer" type="org.hps.recon.ecal.cluster.ReconClusterDriver">
            <logLevel>WARNING</logLevel>
            <outputClusterCollectionName>EcalClusters</outputClusterCollectionName>
        </driver>
        <driver name="CopyCluster" type="org.hps.recon.ecal.cluster.CopyClusterCollectionDriver">
            <inputCollectionName>EcalClusters</inputCollectionName>
            <outputCollectionName>EcalClustersCorr</outputCollectionName>
        </driver>
        <!-- Hodo reconstruction drivers -->
        <driver name="HodoRunningPedestal" type="org.hps.recon.ecal.HodoRunningPedestalDriver">
            <logLevel>CONFIG</logLevel>
        </driver>
        <driver name="HodoRawConverter" type="org.hps.recon.ecal.HodoRawConverterDriver">
           <useRunningPedestal>true</useRunningPedestal>
           <tETAllChannels>8</tETAllChannels>
           <logLevel>CONFIG</logLevel>
        </driver>
        <!-- SVT reconstruction drivers -->
        <!-- Driver used to associate raw tracker hits to corresponding sensor. -->
        <driver name="RawTrackerHitSensorSetup" type="org.lcsim.recon.tracking.digitization.sisim.config.RawTrackerHitSensorSetup">
            <readoutCollections>SVTRawTrackerHits</readoutCollections>
        </driver>
        <!-- Fit the six raw samples and extract the amplitude and t0. --> 
        <driver name="RawTrackerHitFitterDriver" type="org.hps.recon.tracking.RawTrackerHitFitterDriver">
            <fitAlgorithm>Pileup</fitAlgorithm>
            <!--use this to correct for trigger time in MC instead of subtractTriggerTime-->   
            <useTimestamps>false</useTimestamps>     
            <!--offset to get times centered at 0 after timestamp correction-->                 
            <!--<tsCorrectionScale>121</tsCorrectionScale>-->    
            <!--correct for the SVT fit time offset...this should be on if <useTimingConditions> is turned on in readout-->        
            <correctTimeOffset>true</correctTimeOffset>   
            <!--per sensor shift...set false becasue it's not in readout sim-->       
            <correctT0Shift>true</correctT0Shift>      
            <!--use truth time for MC???  typically not used-->             
            <useTruthTime>false</useTruthTime>    
            <!--time of flight corrections-->              
            <subtractTOF>true</subtractTOF>     
            <!--set this false for MC, true for data-->              
            <subtractTriggerTime>true</subtractTriggerTime>           
            <!--per-strip timing correction from database...this should be on i f <useTimingConditions> is turned on in readout  -->  
            <correctChanT0>true</correctChanT0>          
            <isMC>true</isMC>
            <debug>false</debug>
        </driver>
        <!-- 
             Use the fitted raw tracker hits and create 1D clusters using a 
             nearest neighbor algorithm. 
        --> 
        <driver name="TrackerHitDriver" type="org.hps.recon.tracking.DataTrackerHitDriver">
            <neighborDeltaT>50</neighborDeltaT>
        </driver>
        <!-- Create 3D hits from pairs of 1D clusters on adjecent sensors. -->        
        <driver name="HelicalTrackHitDriver" type="org.hps.recon.tracking.HelicalTrackHitDriver">
            <debug>false</debug>
            <clusterTimeCut>100.0</clusterTimeCut>
            <maxDt>100.0</maxDt>
            <clusterAmplitudeCut>400.0</clusterAmplitudeCut>
        </driver>
        <!-- SVT Track finding and fitting. -->
        <driver name="TrackReconSeed123Conf4Extd56" type="org.hps.recon.tracking.TrackerReconDriver">
            <trackCollectionName>Tracks_s123_c4_e56</trackCollectionName>
            <strategyResource>HPS_s123_c4_e56_4hit.xml</strategyResource>
            <debug>false</debug>
            <rmsTimeCut>1000.0</rmsTimeCut>
            <maxTrackerHits>250</maxTrackerHits>
        </driver>
        <driver name="TrackReconSeed123Conf5Extd46" type="org.hps.recon.tracking.TrackerReconDriver">
            <trackCollectionName>Tracks_s123_c5_e46</trackCollectionName>
            <strategyResource>HPS_s123_c5_e46_4hit.xml</strategyResource>
            <debug>false</debug>
            <rmsTimeCut>1000.0</rmsTimeCut>
            <maxTrackerHits>250</maxTrackerHits>
        </driver>
        <driver name="TrackReconSeed567Conf4Extd123" type="org.hps.recon.tracking.TrackerReconDriver">
            <trackCollectionName>Tracks_s567_c4_e123</trackCollectionName>
            <strategyResource>HPS_s567_c4_e123.xml</strategyResource>
            <debug>false</debug>
            <rmsTimeCut>1000.0</rmsTimeCut>
            <maxTrackerHits>250</maxTrackerHits>
        </driver>
        <driver name="TrackReconSeed456Conf3Extd127" type="org.hps.recon.tracking.TrackerReconDriver">
            <trackCollectionName>Tracks_s456_c3_e127</trackCollectionName>
            <strategyResource>HPS_s456_c3_e127.xml</strategyResource>
            <debug>false</debug>
            <rmsTimeCut>1000.0</rmsTimeCut>
            <maxTrackerHits>250</maxTrackerHits>
        </driver>
        <driver name="TrackReconSeed356Conf7Extd124" type="org.hps.recon.tracking.TrackerReconDriver">
            <trackCollectionName>Tracks_s356_c7_e124</trackCollectionName>
            <strategyResource>HPS_s356_c7_e124.xml</strategyResource>
            <debug>false</debug>
            <rmsTimeCut>1000.0</rmsTimeCut>
            <maxTrackerHits>250</maxTrackerHits>
        </driver>
        <driver name="TrackReconSeed235Conf6Extd147" type="org.hps.recon.tracking.TrackerReconDriver">
            <trackCollectionName>Tracks_s235_c6_e147</trackCollectionName>
            <strategyResource>HPS_s235_c6_e147.xml</strategyResource>
            <debug>false</debug>
            <rmsTimeCut>1000.0</rmsTimeCut>
            <maxTrackerHits>250</maxTrackerHits>
        </driver>
        <driver name="TrackReconSeed234Conf6Extd157" type="org.hps.recon.tracking.TrackerReconDriver">
            <trackCollectionName>Tracks_s234_c5_e157</trackCollectionName>
            <strategyResource>HPS_s234_c5_e167_4hit.xml</strategyResource>
            <debug>false</debug>
            <rmsTimeCut>1000.0</rmsTimeCut>
            <maxTrackerHits>250</maxTrackerHits>
        </driver>
        <driver name="TrackClusterTruthMatchingDriver" type="org.hps.recon.utils.TrackClusterTruthMatchingDriver">
            <trackClusterTimeOffset>28</trackClusterTimeOffset>
            <trackCollectionName>KalmanFullTracks</trackCollectionName>
            <trackClusterMatcherAlgo>TrackClusterMatcherMinDistance</trackClusterMatcherAlgo> 
        </driver>

        <driver name="MergeTrackCollections" type="org.hps.recon.tracking.MergeTrackCollections" />
        <driver name="TrackDataDriver" type="org.hps.recon.tracking.TrackDataDriver" />
        <driver name="ReconParticleDriver" type="org.hps.recon.particle.HpsReconParticleDriver" > 
            <ecalClusterCollectionName>EcalClustersCorr</ecalClusterCollectionName>
            <trackCollectionNames>GBLTracks</trackCollectionNames>
            <trackCollectionName>GBLTracks</trackCollectionName>
            <beamPositionX> 0 </beamPositionX>
            <beamSigmaX> 0.05 </beamSigmaX>
            <beamPositionY> 0 </beamPositionY>
            <beamSigmaY> 0.02 </beamSigmaY>
            <beamPositionZ> -7.5 </beamPositionZ>
            <trackClusterTimeOffset>28</trackClusterTimeOffset>
            <maxMatchDt>10</maxMatchDt>
            <useInternalVertexXYPositions>false</useInternalVertexXYPositions>
            <minVertexChisqProb> 0.0 </minVertexChisqProb>
            <maxElectronP>7.0</maxElectronP> 
            <maxVertexP>7.0</maxVertexP>
            <requireClustersForV0>false</requireClustersForV0>
            <trackClusterMatcherAlgo>TrackClusterMatcherMinDistance</trackClusterMatcherAlgo>
            <trackClusterMatchPlots>true</trackClusterMatchPlots>
        </driver>         
        <driver name="ReconParticleDriver_Kalman" type="org.hps.recon.particle.HpsReconParticleDriver" > 
            <ecalClusterCollectionName>EcalClustersCorr</ecalClusterCollectionName>
            <trackCollectionNames>KalmanFullTracks</trackCollectionNames>
            <trackCollectionName>KalmanFullTracks</trackCollectionName>
            <unconstrainedV0CandidatesColName>UnconstrainedV0Candidates_KF</unconstrainedV0CandidatesColName>
            <unconstrainedV0VerticesColName>UnconstrainedV0Vertices_KF</unconstrainedV0VerticesColName>
            <beamConV0CandidatesColName>BeamspotConstrainedV0Candidates_KF</beamConV0CandidatesColName>
            <beamConV0VerticesColName>BeamspotConstrainedV0Vertices_KF</beamConV0VerticesColName>
            <targetConV0CandidatesColName>TargetConstrainedV0Candidates_KF</targetConV0CandidatesColName>
            <targetConV0VerticesColName>TargetConstrainedV0Vertices_KF</targetConV0VerticesColName>
            <finalStateParticlesColName>FinalStateParticles_KF</finalStateParticlesColName>
            <beamPositionX> 0  </beamPositionX>
            <beamSigmaX> 0.05 </beamSigmaX>
            <beamPositionY> 0 </beamPositionY>
            <beamSigmaY> 0.02 </beamSigmaY>
            <beamPositionZ> -7.5 </beamPositionZ>
            <trackClusterTimeOffset>28</trackClusterTimeOffset>
            <maxMatchDt>10</maxMatchDt>
            <useInternalVertexXYPositions>false</useInternalVertexXYPositions>
            <minVertexChisqProb> 0.0 </minVertexChisqProb>
            <maxElectronP>7.0</maxElectronP>
            <maxVertexP>7.0</maxVertexP>
            <requireClustersForV0>false</requireClustersForV0>
            <trackClusterMatcherAlgo>TrackClusterMatcherMinDistance</trackClusterMatcherAlgo>
            <trackClusterMatchPlots>true</trackClusterMatchPlots>
        </driver>         
        <driver name="GBLRefitterDriver" type="org.hps.recon.tracking.gbl.GBLRefitterDriver"/>        
        <driver name="KalmanPatRecDriver" type="org.hps.recon.tracking.kalman.KalmanPatRecDriver">
            <!--<doDebugPlots>false</doDebugPlots>-->
            <verbose>false</verbose>
        </driver>
        <driver name="TrackTruthMatching_KF" type="org.hps.analysis.MC.TrackToMCParticleRelationsDriver">
            <trackCollectionName>KalmanFullTracks</trackCollectionName>
            <kalmanTracks>true</kalmanTracks>
            <debug>false</debug>
        </driver>
        <driver name="TrackTruthMatching_GBL" type="org.hps.analysis.MC.TrackToMCParticleRelationsDriver">
            <trackCollectionName>GBLTracks</trackCollectionName>
            <kalmanTracks>false</kalmanTracks>
            <debug>false</debug>
        </driver>
        <driver name="LCIOWriter" type="org.lcsim.util.loop.LCIODriver">
            <outputFilePath>${outputFile}.slcio</outputFilePath>
        </driver>
        <driver name="CleanupDriver" type="org.lcsim.recon.tracking.digitization.sisim.config.ReadoutCleanupDriver"/>
    </drivers>
</lcsim>
