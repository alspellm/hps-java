//package org.hps.recon.test.ecalScoringPlane;
package org.hps.recon.tracking.kalman;


import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;

import org.hps.recon.ecal.cluster.ClusterUtilities;


import java.util.HashMap;
import java.util.List; 
import java.util.ArrayList; 
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.io.IOException;
import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.TrackData;
import org.hps.util.Pair;
import org.hps.util.RK4integrator;
//import org.lcsim.event.GenericObject;


import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.RawTrackerHit;
//import org.lcsim.event.base.BaseLCRelation;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.util.Driver;
import org.lcsim.geometry.Detector;
import org.lcsim.event.TrackState;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.CalorimeterHit;
//import org.lcsim.event.RawCalorimeterHit;
//import hep.physics.vec.Hep3Vector;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
/** 
 * Driver stolen from Omar to relate a Track to an Ecal scoring plane hit
 *
 **/

public class EcalScoringPlaneDriver extends Driver {

    private org.lcsim.geometry.FieldMap fM;
    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private Map<String, IHistogram2D> plots2D;
    //Histogram special identifiers
    String[] identifiers = {"positive_match", "negative_match", "truth_matched", "no_match"};

    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;
    RelationalTable TrktoData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);

    KFTrackECalClusterMatcher matcher;

    //Counting fake rate
    boolean verbose = true;
    double eleRecoRate;
    double posRecoRate;
    double eleFakeRate=0;
    double posFakeRate=0;
    double eleEfficiency=0;
    double posEfficiency=0;
    double NgeneratedEle = 0;
    double NgeneratedPos = 0;
    double NrecoEle = 0;
    double NrecoPos = 0;
    double NtruthEleClustPairs = 0;
    double NtruthPosClustPairs = 0;
    double eleNegMatch=0;
    double eleTrueMatch=0;
    double posNegMatch=0;
    double posTrueMatch=0;
    double eleNoMatch = 0;
    double posNoMatch = 0;
    double eletotalCount=0;
    double postotalCount=0;
    //time difference for track and cluster
    double trackClusterTimeOffset;

    //Collection Names
    String ecalScoringPlaneHitsCollectionName = "TrackerHitsECal";
    String trackCollectionName = "KalmanFullTracks";
    String trackToScoringPlaneHitRelationsName = "TrackToEcalScoringPlaneHitRelations";

    String ecalClustersCollectionName = "EcalClusters";
    String ecalReadoutHitsCollectionName = "EcalReadoutHits";
    String ecalTruthRelationsName = "EcalTruthRelations";

    private Set<SimTrackerHit> simhitsontrack = new HashSet<SimTrackerHit>();


    public void setTrackClusterTimeOffset(double input) {
        trackClusterTimeOffset = input;
    }
    public void setTrackCollectionName(String trackCollectionName) {
        this.trackCollectionName = trackCollectionName;
    }

    public void saveHistograms() {
        System.out.println("[EcalScoringPlaneDriver] Saving Histogram for " + this.trackCollectionName);
        String rootFile = String.format("%s_truthTrackClusterMatching.root",this.trackCollectionName);
        RootFileStore store = new RootFileStore(rootFile);
        try {
            store.open();
            store.add(tree);
            store.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void bookHistograms() {

        String trackType = this.trackCollectionName;
        System.out.println("[EcalScoringPlaneDriver] Booking Histograms for " + trackType);
        plots1D = new HashMap<String, IHistogram1D>();
        plots2D = new HashMap<String, IHistogram2D>();
        tree = IAnalysisFactory.create().createTreeFactory().create();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

        //truth track position at scoring plane
        plots1D.put(String.format("%s_ele_truth_extrapToScoringPlane_z",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_extrapToScoringPlane_z",this.trackCollectionName), 1500, 0, 1500));
        plots1D.put(String.format("%s_pos_truth_extrapToScoringPlane_z",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_extrapToScoringPlane_z",this.trackCollectionName), 1500, 0, 1500));
        //Checking the Z position of TrackerHitsEcal collection. Should all be
        //at same Z position
        plots1D.put(String.format("%s_TrackerHitsEcal_z",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_TrackerHitsEcal_z",this.trackCollectionName), 1500, 0, 1500));

        plots1D.put(String.format("%s_ele_Track_atEcal_ScoringPlane_dx",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Track_atEcal_ScoringPlane_dx",this.trackCollectionName), 200, -200, 200));
        plots1D.put(String.format("%s_ele_Track_atEcal_ScoringPlane_dy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Track_atEcal_ScoringPlane_dy",this.trackCollectionName), 200, -200, 200));
        plots1D.put(String.format("%s_ele_Track_atEcal_ScoringPlane_dz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Track_atEcal_ScoringPlane_dz",this.trackCollectionName), 200, -200, 200));
        plots1D.put(String.format("%s_ele_Track_atEcal_ScoringPlane_dr",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Track_atEcal_ScoringPlane_dr",this.trackCollectionName), 200, -200, 200));
        plots1D.put(String.format("%s_ele_Track_atEcal_ScoringPlane_dt",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Track_atEcal_ScoringPlane_dt",this.trackCollectionName), 200, -200, 200));

        plots1D.put(String.format("%s_pos_Track_atEcal_ScoringPlane_dx",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Track_atEcal_ScoringPlane_dx",this.trackCollectionName), 200, -200, 200));
        plots1D.put(String.format("%s_pos_Track_atEcal_ScoringPlane_dy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Track_atEcal_ScoringPlane_dy",this.trackCollectionName), 200, -200, 200));
        plots1D.put(String.format("%s_pos_Track_atEcal_ScoringPlane_dz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Track_atEcal_ScoringPlane_dz",this.trackCollectionName), 200, -200, 200));
        plots1D.put(String.format("%s_pos_Track_atEcal_ScoringPlane_dr",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Track_atEcal_ScoringPlane_dr",this.trackCollectionName), 200, -200, 200));
        plots1D.put(String.format("%s_pos_Track_atEcal_ScoringPlane_dt",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Track_atEcal_ScoringPlane_dt",this.trackCollectionName), 200, -200, 200));
        //Track extrapolation to Ecal: Momentum vs truth-extrap position
        //residuals
        plots2D.put(String.format("%s_ele_track_truth-extrapolation_dx_v_P",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_track_truth-extrapolation_dx_v_P",this.trackCollectionName), 600, -300, 300,300,0,3));
        plots2D.put(String.format("%s_ele_track_truth-extrapolation_dy_v_P",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_track_truth-extrapolation_dy_v_P",this.trackCollectionName), 600, -300, 300,300,0,3));
        plots2D.put(String.format("%s_ele_truthMCP_RK4ScoringPlaneToEcal_ZvP",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_truthMCP_RK4ScoringPlaneToEcal_ZvP",this.trackCollectionName), 1500, 0, 1500,300,0,3));

        plots2D.put(String.format("%s_pos_track_truth-extrapolation_dx_v_P",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_track_truth-extrapolation_dx_v_P",this.trackCollectionName), 600, -300, 300,300,0,3));
        plots2D.put(String.format("%s_pos_track_truth-extrapolation_dy_v_P",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_track_truth-extrapolation_dy_v_P",this.trackCollectionName), 600, -300, 300,300,0,3));
        plots2D.put(String.format("%s_pos_truthMCP_RK4ScoringPlaneToEcal_ZvP",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_truthMCP_RK4ScoringPlaneToEcal_ZvP",this.trackCollectionName), 1500, 0, 1500,300,0,3));


        for(String id: identifiers){
            plots1D.put(String.format("%s_ele_track_cluster_%s_dx",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_ele_track_cluster_%s_dx",this.trackCollectionName,id), 200, -200, 200));
            plots1D.put(String.format("%s_ele_track_cluster_%s_dy",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_ele_track_cluster_%s_dy",this.trackCollectionName,id), 200, -200, 200));
            plots1D.put(String.format("%s_ele_track_cluster_%s_dz",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_ele_track_cluster_%s_dz",this.trackCollectionName,id), 200, -200, 200));
            plots1D.put(String.format("%s_ele_track_cluster_%s_dr",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_ele_track_cluster_%s_dr",this.trackCollectionName,id), 200, -200, 200));
            plots1D.put(String.format("%s_ele_track_cluster_%s_dt",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_ele_track_cluster_%s_dt",this.trackCollectionName,id), 200, -200, 200));

            plots1D.put(String.format("%s_pos_track_cluster_%s_dx",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_pos_track_cluster_%s_dx",this.trackCollectionName,id), 200, -200, 200));
            plots1D.put(String.format("%s_pos_track_cluster_%s_dy",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_pos_track_cluster_%s_dy",this.trackCollectionName,id), 200, -200, 200));
            plots1D.put(String.format("%s_pos_track_cluster_%s_dz",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_pos_track_cluster_%s_dz",this.trackCollectionName,id), 200, -200, 200));
            plots1D.put(String.format("%s_pos_track_cluster_%s_dr",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_pos_track_cluster_%s_dr",this.trackCollectionName,id), 200, -200, 200));
            plots1D.put(String.format("%s_pos_track_cluster_%s_dt",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_pos_track_cluster_%s_dt",this.trackCollectionName,id), 200, -200, 200));


            //plots1D.put(String.format("%s_cluster_energy_%s",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_cluster_energy_%s",this.trackCollectionName,id),1000,-1000,1000));
            //plots1D.put(String.format("%s_track_momentum_%s",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_track_momentum_%s",this.trackCollectionName,id),1000,-1000,1000));
            //plots1D.put(String.format("%s_momentum_energy_residual_%s",this.trackCollectionName,id), histogramFactory.createHistogram1D(String.format("%s_momentum_energy_residual%s",this.trackCollectionName,id),1000,-1000,1000));
        }

        //Matched vs Unmatched Track positions   
        /*
        plots1D.put(String.format("%s_ele_track_NoMatchedCluster_x",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_NoMatchedCluster_x",this.trackCollectionName), 1000, -1000, 1000));
        plots1D.put(String.format("%s_ele_track_NoMatchedCluster_y",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_NoMatchedCluster_y",this.trackCollectionName), 1000, -1000, 1000));
        plots1D.put(String.format("%s_ele_track_NoMatchedCluster_z",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_NoMatchedCluster_z",this.trackCollectionName), 2000, -2000, 2000));
        plots1D.put(String.format("%s_ele_track_NoMatchedCluster_t",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_NoMatchedCluster_t",this.trackCollectionName), 2000, -2000, 2000));
        plots1D.put(String.format("%s_pos_track_NoMatchedCluster_x",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_NoMatchedCluster_x",this.trackCollectionName), 1000, -1000, 1000));
        plots1D.put(String.format("%s_pos_track_NoMatchedCluster_y",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_NoMatchedCluster_y",this.trackCollectionName), 1000, -1000, 1000));
        plots1D.put(String.format("%s_pos_track_NoMatchedCluster_z",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_NoMatchedCluster_z",this.trackCollectionName), 2000, -2000, 2000));
        plots1D.put(String.format("%s_pos_track_NoMatchedCluster_t",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_NoMatchedCluster_t",this.trackCollectionName), 2000, -2000, 2000));
        */


        /*
        plots1D.put(String.format("%s_ele_track_MatchedCluster_x",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_MatchedCluster_x",this.trackCollectionName), 1000, -1000, 1000));
        plots1D.put(String.format("%s_ele_track_MatchedCluster_y",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_MatchedCluster_y",this.trackCollectionName), 1000, -1000, 1000));
        plots1D.put(String.format("%s_ele_track_MatchedCluster_z",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_MatchedCluster_z",this.trackCollectionName), 2000, -2000, 2000));
        plots1D.put(String.format("%s_ele_track_MatchedCluster_t",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_MatchedCluster_t",this.trackCollectionName), 2000, -2000, 2000));
        plots1D.put(String.format("%s_pos_track_MatchedCluster_x",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_MatchedCluster_x",this.trackCollectionName), 1000, -1000, 1000));
        plots1D.put(String.format("%s_pos_track_MatchedCluster_y",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_MatchedCluster_y",this.trackCollectionName), 1000, -1000, 1000));
        plots1D.put(String.format("%s_pos_track_MatchedCluster_z",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_MatchedCluster_z",this.trackCollectionName), 2000, -2000, 2000));
        plots1D.put(String.format("%s_pos_track_MatchedCluster_t",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_MatchedCluster_t",this.trackCollectionName), 2000, -2000, 2000));
        */
        ///
        plots1D.put(String.format("%s_ele_fakeRate",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_fakeRate",this.trackCollectionName), 10, -1, 10));
        plots1D.put(String.format("%s_pos_fakeRate",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_fakeRate",this.trackCollectionName), 10, -1, 10));
        plots1D.put(String.format("%s_ele_Efficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Efficiency",this.trackCollectionName), 10, -1, 10));
        plots1D.put(String.format("%s_pos_Efficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Efficiency",this.trackCollectionName), 10, -1, 10));
        plots1D.put(String.format("%s_ele_Reco_Rate",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Reco_Rate",this.trackCollectionName), 10, -1, 10));
        plots1D.put(String.format("%s_pos_Reco_Rate",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Reco_Rate",this.trackCollectionName), 10, -1, 10));


        //Hit multiplicity for truth matching
        plots1D.put(String.format("%s_Track_max_mcParticleMultiplicity",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_Track_max_mcParticleMultiplicity",this.trackCollectionName),30, 0, 30));

        //2D Plots
        plots2D.put(String.format("%s_ele_Track_Pos_at_Ecal; x mm; y mm",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_Track_Pos_at_Ecal; x mm; y mm",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_Track_Pos_at_Ecal; x mm; y mm",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_Track_Pos_at_Ecal; x mm; y mm",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        //Cluster Positions
        //plots2D.put(String.format("ele_EcalCluster_position; x mm; y mm",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("ele_EcalCluster_position; x mm; y mm",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("EcalCluster_position; x mm; y mm",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("EcalCluster_position; x mm; y mm",this.trackCollectionName),1000, -500, 500,1000, -500, 500));



    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void startOfData() {
        System.out.println("Starting job");
        bookHistograms();
        matcher = new KFTrackECalClusterMatcher(this.trackCollectionName);
    }

    public void endOfData() {
        eleFakeRate = eleNegMatch/(eleNegMatch+eleTrueMatch);
        posFakeRate = posNegMatch/(posNegMatch + posTrueMatch);
        eleEfficiency = (eleNegMatch+eleTrueMatch)/NtruthEleClustPairs;
        posEfficiency = (posNegMatch+posTrueMatch)/NtruthPosClustPairs;
        eleRecoRate = NrecoEle/NgeneratedEle;
        posRecoRate = NrecoPos/NgeneratedPos;

        for(int i=0; i < Math.round(eleFakeRate*100); i++)
            plots1D.get(String.format("%s_ele_fakeRate",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(posFakeRate*100); i++)
            plots1D.get(String.format("%s_pos_fakeRate",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(posEfficiency*100); i++)
            plots1D.get(String.format("%s_pos_Efficiency",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(eleEfficiency*100); i++)
            plots1D.get(String.format("%s_ele_Efficiency",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(eleRecoRate*100); i++)
            plots1D.get(String.format("%s_ele_Reco_Rate",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(posRecoRate*100); i++)
            plots1D.get(String.format("%s_ele_Reco_Rate",this.trackCollectionName)).fill(1.0);
        saveHistograms();
    }

    protected void process(EventHeader event) {



        //Get collection of Ecal clusters
        //If event has no collection of tracks, skip
        if(!event.hasCollection(Track.class, trackCollectionName)) return;
        //If even doesnt have collection of Ecal scoring plane hits, skip
        if(!event.hasCollection(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName)) return;

        //Get EcalClusters from event
        List<Cluster> clusters = event.get(Cluster.class, ecalClustersCollectionName);
        for(Cluster cluster : clusters) {
            double clustx = cluster.getPosition()[0];
            double clusty = cluster.getPosition()[1];
            double clustz = cluster.getPosition()[2];
            plots2D.get(String.format("EcalCluster_position; x mm; y mm")).fill(clustx,clusty);
        }
        // Get collection of tracks from event
        List<Track> tracks = event.get(Track.class, trackCollectionName);
        List<Track> truthmatchedTracks = new ArrayList<Track>();
        Map<Track, Cluster> truthTrackClusterMap = new HashMap<Track, Cluster>();

        hitToRotated = TrackUtils.getHitToRotatedTable(event);
        hitToStrips = TrackUtils.getHitToStripsTable(event);
        List<TrackData> TrackData;
        List<LCRelation> trackRelations;
        TrackData trackdata;
        if (this.trackCollectionName.contains("KalmanFullTracks")) {
            TrackData = event.get(TrackData.class, "KFTrackData");
            trackRelations = event.get(LCRelation.class, "KFTrackDataRelations");
            for (LCRelation relation : trackRelations) {
                if (relation != null && relation.getTo() != null){
                    TrktoData.add(relation.getFrom(), relation.getTo());
                }
            }
        }

        System.out.println("[EcalScoringPlane] track length: " + tracks.size()); 
        for(Track track : tracks){
            SimTrackerHit matchedScoringplaneHit = getTrackScoringPlaneHit(event, track, ecalScoringPlaneHitsCollectionName);
            if(matchedScoringplaneHit != null)
                this.trackScoringPlanePlots(event, track, matchedScoringplaneHit);

            Cluster truthmatchedCluster = truthMatchTracktoCluster(event, track, clusters); 
            if(truthmatchedCluster != null){
                truthmatchedTracks.add(track);
                truthTrackClusterMap.put(track, truthmatchedCluster);
            }
        }

        //Use matching algorithm to match tracks to clusters (without truth
        //info)
        Map<Track, Cluster> matchedTrackClusterMap = new HashMap<Track,Cluster>();
        matchedTrackClusterMap = matcher.trackClusterMatcher(truthmatchedTracks,event, this.trackCollectionName, clusters, this.trackClusterTimeOffset);

        System.out.println("[EcalScoringPlane] truthmatched track length: " + truthmatchedTracks.size()); 
        if(truthmatchedTracks != null){
            for(Track track : truthmatchedTracks) {
                int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
                //Check purity and fake rate of KFTrackECalClusterMatcher
                if(charge < 0)
                    NtruthEleClustPairs = NtruthEleClustPairs + 1;
                else
                    NtruthPosClustPairs = NtruthPosClustPairs + 1;

                checkTrackClusterMatch(track, matchedTrackClusterMap.get(track), truthTrackClusterMap.get(track));
                trackClusterAnalysis(track,truthTrackClusterMap.get(track),"truth_matched");
            }
        }

    }

    public Cluster truthMatchTracktoCluster(EventHeader event, Track track, List<Cluster> clusters){

        /**
         * Get the MC particle associated with a track.
         * Fill mape with Track -> MCParticle
        **/
        String MCHitInputCollectionName = "TrackerHits";
        List<SimTrackerHit> allsimhits = event.get(SimTrackerHit.class, MCHitInputCollectionName);
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, "SVTTrueHitRelations")) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, "SVTTrueHitRelations");
            for (LCRelation relation : trueHitRelations)
                if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                    rawtomc.add(relation.getFrom(), relation.getTo());
        }

        boolean isKalman = false;
        if(this.trackCollectionName.contains("KalmanFullTracks"))
            isKalman = true;
  
        //create instance of TrackTruthMatching (by Matt Solt)
        TrackTruthMatching truthmatcher = new TrackTruthMatching(track, rawtomc, allsimhits, isKalman);
        MCParticle trackMCP = truthmatcher.getMCParticle();
        //MCParticle trackMCP = this.getTrackMCP(track,event);

        if(trackMCP == null) return null;

        /**
         * Truth matching Tracks with Clusters via track_MCParticle ==
         * cluster_MCParticle.
         * */

        Cluster trackMatchedCluster = null;
        for (Cluster cluster : clusters) {
            MCParticle clusterMCP = getClusterMCP(cluster, event);
            //if MCParticles of cluster and track match, map the two.
            if(clusterMCP == trackMCP) {
                trackMatchedCluster = cluster;                    
                break;
            }
        }
        return trackMatchedCluster;
    }

    public SimTrackerHit getTrackScoringPlaneHit(EventHeader event, Track track, String ecalScoringPlaneHitsCollectionName) {

        List<SimTrackerHit> scoringPlaneHits = event.get(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName);

        MCParticle trackMCP = this.getTrackMCP(track,event);

        SimTrackerHit matchedScoringPlaneHit = null;
        for(SimTrackerHit scoringPlaneHit : scoringPlaneHits){
            // If the MC particles don't match, move on to the next particle
            if(!(scoringPlaneHit.getMCParticle() == trackMCP)) continue;
            matchedScoringPlaneHit = scoringPlaneHit;
            // Once a match is found, there is no need to loop through the rest of the list
            break;
        }
        return matchedScoringPlaneHit;
    }

    public void checkTrackClusterMatch(Track track,Cluster matchedCluster, Cluster truthmatchedCluster){
        
        //Get track parameters
        int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());

        if(matchedCluster == null){
            trackClusterAnalysis(track,matchedCluster,"no_match");
            if(charge > 0) 
                posNoMatch = posNoMatch + 1.0;
            else
                eleNoMatch = eleNoMatch + 1.0;
            return;
        }

        if(matchedCluster == truthmatchedCluster){
            trackClusterAnalysis(track, matchedCluster,"positive_match");
            if(charge > 0)
                posTrueMatch = posTrueMatch + 1.0;
            else
                eleTrueMatch = eleTrueMatch + 1.0;
        }
        else{
            trackClusterAnalysis(track, matchedCluster,"negative_match");
            if(charge > 0) 
                posNegMatch = posNegMatch + 1.0;
            else
                eleNegMatch = eleNegMatch + 1.0;
        }

        eletotalCount = eleTrueMatch + eleNegMatch + eleNoMatch;
        postotalCount = posTrueMatch + posNegMatch + posNoMatch;
    }

    

    public void trackScoringPlanePlots(EventHeader event, Track track, SimTrackerHit scoringplaneHit) {

        //Comparison plots for Track extrapolated to ECal vs Truth ScoringPlaneHit at Ecal
    
        double trackx;
        double tracky;
        double trackz;
        double truthxpos;
        double truthypos;
        double truthzpos;
        double dxoffset;

        double trackT;
        double simTrackT;

        if (this.trackCollectionName.contains("GBLTracks")){
            trackT = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
        }
        else {
            TrackData trackdata = (TrackData) TrktoData.from(track);
            trackT = trackdata.getTrackTime();
        }
        simTrackT = scoringplaneHit.getTime();

        //Make histograms of truth vs extrapolation
        if(this.trackCollectionName.contains("GBLTracks")) {
            trackx = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[1];
            tracky = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[2];
            trackz = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[0];
            //dxoffset = -4.1;
            dxoffset = 0.0;
        }

        else {
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trackx = ts_ecalPos[0];
            tracky = ts_ecalPos[1];
            trackz = ts_ecalPos[2];
            dxoffset = 0.0;
        }

        truthxpos = scoringplaneHit.getPoint()[0];
        truthypos = scoringplaneHit.getPoint()[1];
        truthzpos = scoringplaneHit.getPoint()[2];

        //multiply charge by factor of -1 (WHY!?)
        int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());

        double[] trackP = track.getMomentum();
        double trackPmag = Math.sqrt(Math.pow(trackP[0],2) + Math.pow(trackP[1],2) + Math.pow(trackP[2],2));
        double[] truthP = scoringplaneHit.getMomentum();


        //In PhysicsRun2016 detector geometry, scoring plane is located
        //upstream of Ecal face by ~30 mm...
        //Need to extrapolate MCP scoring plane hits to the Ecal, so that we
        //can compare track positions at Ecal with MCP hit at Ecal
        //rotate HPS coordinates to Tracking coordinates (XYZ)HPS -> (ZXY)
        Hep3Vector truthPVecTracking = new BasicHep3Vector(truthP[2],truthP[0],truthP[1]);
        Hep3Vector truthPosVecTracking = new BasicHep3Vector(truthzpos,truthxpos,truthypos);

        Pair<Hep3Vector,Hep3Vector> extrapVecsTracking = null;
        Detector detector = event.getDetector();
        fM = detector.getFieldMap();
        RK4integrator RKint = new RK4integrator(charge,1, fM);
        extrapVecsTracking = RKint.integrate(truthPosVecTracking,truthPVecTracking,1393-scoringplaneHit.getPoint()[2]);
        Hep3Vector RKextrapPos = extrapVecsTracking.getFirstElement();
        truthxpos = RKextrapPos.y();
        truthypos = RKextrapPos.z();
        truthzpos = RKextrapPos.x();

        //Take position residuals for track and extrap truth hits at ecal
        double dx = truthxpos - trackx;
        double dy = truthypos - tracky;
        double dz = truthzpos - trackz;
        double dr = Math.sqrt(Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2));
        double dt = simTrackT - trackT;

        //Make plots
        if(charge < 0) {


            //Track X,Y position plots at Ecal
            plots2D.get(String.format("%s_ele_Track_Pos_at_Ecal; x mm; y mm",this.trackCollectionName)).fill(trackx,tracky);
            plots1D.get(String.format("%s_ele_truth_extrapToScoringPlane_z",this.trackCollectionName)).fill(truthzpos);
            //Extrapolated Track momentum vs truth position residuals
            plots2D.get(String.format("%s_ele_track_truth-extrapolation_dx_v_P",this.trackCollectionName)).fill(dx, trackPmag);
            plots2D.get(String.format("%s_ele_track_truth-extrapolation_dy_v_P",this.trackCollectionName)).fill(dy, trackPmag);
            plots2D.get(String.format("%s_ele_truthMCP_RK4ScoringPlaneToEcal_ZvP",this.trackCollectionName)).fill(truthzpos, trackPmag);

            //Track at Ecal vs Cluster residuals
            plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dx",this.trackCollectionName)).fill(dx);
            plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dy",this.trackCollectionName)).fill(dy);
            plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dz",this.trackCollectionName)).fill(dz);
            plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dr",this.trackCollectionName)).fill(dr);
            plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dt",this.trackCollectionName)).fill(dt);
        }
        else {
            plots1D.get(String.format("%s_pos_truth_extrapToScoringPlane_z",this.trackCollectionName)).fill(truthzpos);
            //Track X,Y position at Ecal
            plots2D.get(String.format("%s_pos_Track_Pos_at_Ecal; x mm; y mm",this.trackCollectionName)).fill(trackx,tracky);
            //Extrapolated Track P vs truth position residuals
            plots2D.get(String.format("%s_pos_track_truth-extrapolation_dx_v_P",this.trackCollectionName)).fill(dx, trackPmag);
            plots2D.get(String.format("%s_pos_track_truth-extrapolation_dy_v_P",this.trackCollectionName)).fill(dy, trackPmag);
            plots2D.get(String.format("%s_pos_truthMCP_RK4ScoringPlaneToEcal_ZvP",this.trackCollectionName)).fill(truthzpos, trackPmag);

            //Track vs Cluster residuals at Ecal
            plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dx",this.trackCollectionName)).fill(dx);
            plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dy",this.trackCollectionName)).fill(dy);
            plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dz",this.trackCollectionName)).fill(dz);
            plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dr",this.trackCollectionName)).fill(dr);
            plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dt",this.trackCollectionName)).fill(dt);
        }

    }

    public void trackClusterAnalysis(Track track, Cluster cluster,  String identifier) {
        //For comparing extrapolated track position  with position of
        //truth-matched cluster
        int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
        double trackx;
        double tracky;
        double trackz;
        double clustx;
        double clusty;
        double clustz;
        double dxoffset;
        double clusTime;
        double dt;
        double dx;
        double dy;
        double dz;
        double dr;
        double clusterEnergy; 
        String id = identifier; //positive_match, negative_match, truth_match
        double trackT;
        if (this.trackCollectionName.contains("GBLTracks")){
            trackT = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
            trackx = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[1];
            tracky = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[2];
            trackz = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[0];
          //  dxoffset = -5.5;
            dxoffset = 0.0;
        }
        else {
            TrackData trackdata = (TrackData) TrktoData.from(track);
            trackT = trackdata.getTrackTime();
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trackx = ts_ecalPos[0];
            tracky = ts_ecalPos[1];
            trackz = ts_ecalPos[2];
            dxoffset = 0.0;
        }
        //double tracktOffset = 4.0;
        //double[] trackP = new double[3];
        //double[] trackP = track.getMomentum();
        
        if(cluster == null){
            clustx = 0;
            clusty = 0;
            clustz = 0;
            clusTime = 0;
            clusterEnergy = 0.0;
            dt = trackT;
        }
        else {
            clustx = cluster.getPosition()[0];
            clusty = cluster.getPosition()[1];
            clustz = cluster.getPosition()[2];
            clusterEnergy = cluster.getEnergy();
            clusTime = ClusterUtilities.getSeedHitTime(cluster);
            dt = clusTime - trackClusterTimeOffset - trackT;
        }

        dx = trackx - clustx;
        dy = tracky - clusty;
        dz = trackz - clustz;
        dr = Math.sqrt(Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2));

        //Make plots
        if(charge < 0) {
            if(track == null)
            plots1D.get(String.format("%s_ele_track_cluster_%s_dx",this.trackCollectionName,id)).fill(dx);
            plots1D.get(String.format("%s_ele_track_cluster_%s_dy",this.trackCollectionName,id)).fill(dy);
            plots1D.get(String.format("%s_ele_track_cluster_%s_dz",this.trackCollectionName,id)).fill(dz);
            plots1D.get(String.format("%s_ele_track_cluster_%s_dr",this.trackCollectionName,id)).fill(dr);
            plots1D.get(String.format("%s_ele_track_cluster_%s_dt",this.trackCollectionName,id)).fill(dt);
        }
        else {
            plots1D.get(String.format("%s_pos_track_cluster_%s_dx",this.trackCollectionName,id)).fill(dx);
            plots1D.get(String.format("%s_pos_track_cluster_%s_dy",this.trackCollectionName,id)).fill(dy);
            plots1D.get(String.format("%s_pos_track_cluster_%s_dz",this.trackCollectionName,id)).fill(dz);
            plots1D.get(String.format("%s_pos_track_cluster_%s_dr",this.trackCollectionName,id)).fill(dr);
            plots1D.get(String.format("%s_pos_track_cluster_%s_dt",this.trackCollectionName,id)).fill(dt);
        }

    }

/**
     * Get the MC particle associated with a track.
     * 
     * @param track : Track to get the MC particle for
     * @return The MC particle associated with the track
     */
    private MCParticle getClusterMCP(Cluster cluster, EventHeader event) {
        List <RawTrackerHit> ecalReadoutHits = event.get(RawTrackerHit.class,ecalReadoutHitsCollectionName);
        Map <MCParticle, int[]>mcParticleMultiplicity = new HashMap<MCParticle, int[]>();
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, ecalTruthRelationsName )) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, ecalTruthRelationsName);
            for (LCRelation relation : trueHitRelations)
                if (relation != null && relation.getFrom() != null && relation.getTo() != null){
                    rawtomc.add(relation.getFrom(), relation.getTo());
                }
                else System.out.println("failed to build cluster mc relation");
        }

        //Match Cluster->CalorimeterHit(s) to corresponding EcalReadoutHits via
        //CellID match
        List<CalorimeterHit> calorimeterhits = cluster.getCalorimeterHits();
        CalorimeterHit seedhit = ClusterUtilities.findSeedHit(cluster);
        RawTrackerHit readoutMatchHit = null;
        long cellID = seedhit.getCellID();
        for(RawTrackerHit ecalReadoutHit : ecalReadoutHits) {
            long recellID = ecalReadoutHit.getCellID();
            if(cellID == recellID) {
                readoutMatchHit = ecalReadoutHit;
            }
        }

        Set<SimCalorimeterHit> simcalhits = rawtomc.allFrom(readoutMatchHit);
        double maxMCPEnergy = 0.0;
        MCParticle largestEnergyMCP = null;
        for(SimCalorimeterHit simcalhit : simcalhits) {
            for(int i=0; i < simcalhit.getMCParticleCount(); i++){
                if(simcalhit.getMCParticle(i).getEnergy() > maxMCPEnergy) {
                    maxMCPEnergy = simcalhit.getMCParticle(i).getEnergy();
                    largestEnergyMCP = simcalhit.getMCParticle(i);
                }
            }
        }

        return largestEnergyMCP;
    }
        

    private MCParticle getTrackMCP(Track track, EventHeader event){
        
        Map <MCParticle, int[]>mcParticleMultiplicity = new HashMap<MCParticle, int[]>();

        //Retrieve rawhits to mc
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, "SVTTrueHitRelations")) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, "SVTTrueHitRelations");
            for (LCRelation relation : trueHitRelations)
                if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                    rawtomc.add(relation.getFrom(), relation.getTo());
        }

        //Retrieve all simulated hits
        String MCHitInputCollectionName = "TrackerHits";
        List<SimTrackerHit> allsimhits = event.get(SimTrackerHit.class, MCHitInputCollectionName);
        MCParticle particle;
        for(TrackerHit hit : track.getTrackerHits()){
            int rawHitsLength = hit.getRawHits().size();
            List<RawTrackerHit> rawhits = hit.getRawHits();
            SimTrackerHit largestSimhit = null;
            for(RawTrackerHit rawhit : rawhits){
                double simhitMaxE = 0.0;
                Set<SimTrackerHit> simhits = rawtomc.allFrom(rawhit);
                for(SimTrackerHit simhit : simhits){
                    if (simhit != null && simhit.getMCParticle() != null) {
                        double simhitEnergy = simhit.getdEdx();
                        if(simhitEnergy > simhitMaxE){
                            simhitMaxE = simhitEnergy;
                            largestSimhit = simhit;
                        }
                    }
                }
                            

                /*IDDecoder decoder = rawhit.getIDDecoder();
                decoder.setID(rawhit.getCellID());
                int layer = decoder.getLayer();
                SymmetricMatrix covmatrix = new SymmetricMatrix(3,hit.getCovMatrix(),false);
                BarrelEndcapFlag bef = BarrelEndcapFlag.UNKNOWN;
                BasicHep3Vector gpos = new BasicHep3Vector(hit.getPosition());
                */
                SimTrackerHit simhit = largestSimhit;
                if (simhit != null && simhit.getMCParticle() != null) {
                    simhitsontrack.add(simhit);
                    particle = simhit.getMCParticle();
                    if(!mcParticleMultiplicity.containsKey(particle)){
                        mcParticleMultiplicity.put(particle, new int[1]);
                        mcParticleMultiplicity.get(particle)[0] = 0;
                    }

                    mcParticleMultiplicity.get(particle)[0]++;

                    break;
                }
            }

        }

        // Look for the MC particle that occurs the most of the track
        int maxValue = 0;
        particle = null;
        for(Map.Entry<MCParticle, int[]> entry : mcParticleMultiplicity.entrySet()){
            if(maxValue < entry.getValue()[0]){
                particle = entry.getKey();
                maxValue = entry.getValue()[0];
            }
        }
        plots1D.get(String.format("%s_Track_max_mcParticleMultiplicity",this.trackCollectionName)).fill(maxValue);
        if(this.trackCollectionName.contains("GBLTracks")){
            if(maxValue > 4)
                return particle;
            else
                return null;
        }
        else if(this.trackCollectionName.contains("KalmanFullTracks")){
            if(maxValue > 9)
                return particle;
            else
                return null;
        }
        else
            return null;
    }

    public void getTrackRecoEfficiency(EventHeader event){
        /*
        Map <MCParticle, int[]>mcParticleMultiplicity = new HashMap<MCParticle, int[]>();
        List<SimTrackerHit> simhits = event.get(SimTrackerHit.class,MCHitInputCollectionName);
        List<MCParticle> mcParticles = new ArrayList<MCParticle>();
        List<MCParticle> generatedEles = new ArrayList<MCParticle>();
        List<MCParticle> generatedPoss = new ArrayList<MCParticle>();
        for(SimTrackerHit simhit : simhits){
            if(simhit.getMCParticle == null)
                continue;
           // mcParticles.add(simhit.getMCParticle());
            MCParticle particle = simhit.getMCParticle();
            if(!mcParticleMultiplicity.containsKey(particle)){
                mcParticleMultiplicity.put(particle, new int[1]);
                mcParticleMultiplicity.get(particle)[0] = 0;
            }
            mcParticleMultiplicity.get(particle)[0]++;


        }

        for(Map.Entry<MCParticle, int[]> entry : mcParticleMultiplicity.entrySet()){
            particle = entry.getKey();
            double hitcount = entry.getValue()[0];
            int charge =  particle.getCharge();
            if(this.trackCollectionName.contains("GBLTracks")){
                if(hitcount > 4){
                    if(charge < 0) 
                        generatedEles.add(particle);
                    else
                        generatedPoss.add(particle);
                }
            }
            else if(this.trackCollectionName.contains("KalmanFullTracks")){
                if(hitcount > 9){
                    if(charge < 0)
                        generatedEles.add(particle);
                    else
                        generatedPoss.add(particle);
                }
            }

        }
        */
    }
}


