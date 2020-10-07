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
//import org.lcsim.event.GenericObject;


import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.base.BaseLCRelation;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.util.Driver;
import org.lcsim.geometry.IDDecoder;
import org.lcsim.geometry.subdetector.BarrelEndcapFlag;
import hep.physics.vec.BasicHep3Vector;
import org.lcsim.event.TrackState;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.CalorimeterHit;
//import org.lcsim.event.RawCalorimeterHit;
//import hep.physics.vec.Hep3Vector;
import hep.physics.matrix.SymmetricMatrix;
/** 
 * Driver stolen from Omar to relate a Track to an Ecal scoring plane hit
 *
 **/

public class EcalScoringPlaneDriver extends Driver {

    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private Map<String, IHistogram2D> plots2D;
    //Histogram special identifiers
    String[] identifiers = {"positive_match", "negative_match", "truth_matched"};

    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;
    KFTrackECalClusterMatcher matcher;

    //Counting fake rate
    boolean verbose = true;
    double eleRecoRate;
    double posRecoRate;
    double eleFakeRate=0;
    double posFakeRate=0;
    double eleInefficiency=0;
    double posInefficiency=0;
    double NgeneratedEle = 0;
    double NgeneratedPos = 0;
    double NrecoEle = 0;
    double NrecoPos = 0;
    double NrecoElewMCP = 0;
    double NrecoPoswMCP = 0;
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
    String trackToMCParticleRelationsName = "TrackToMCParticleRelations";

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
        plots1D = new HashMap<String, IHistogram1D>();
        plots2D = new HashMap<String, IHistogram2D>();
        tree = IAnalysisFactory.create().createTreeFactory().create();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

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
        plots1D.put(String.format("%s_ele_Inefficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Inefficiency",this.trackCollectionName), 10, -1, 10));
        plots1D.put(String.format("%s_pos_Inefficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Inefficiency",this.trackCollectionName), 10, -1, 10));
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
        eleFakeRate = eleNegMatch/NrecoElewMCP;
        posFakeRate = posNegMatch/NrecoPoswMCP;
        eleInefficiency = eleNoMatch/NrecoElewMCP;
        posInefficiency = posNoMatch/NrecoPoswMCP;
        eleRecoRate = NrecoEle/NgeneratedEle;
        posRecoRate = NrecoPos/NgeneratedPos;

        for(int i=0; i < Math.round(eleFakeRate*100); i++)
            plots1D.get(String.format("%s_ele_fakeRate",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(posFakeRate*100); i++)
            plots1D.get(String.format("%s_pos_fakeRate",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(posInefficiency*100); i++)
            plots1D.get(String.format("%s_pos_Inefficiency",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(eleInefficiency*100); i++)
            plots1D.get(String.format("%s_ele_Inefficiency",this.trackCollectionName)).fill(1.0);
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
        //Get collection of Ecal scoring plane hits from event
        List<SimTrackerHit> scoringPlaneHits = event.get(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName);

        //Create collection to hold scoring plane hits found to match a track
        List<SimTrackerHit> matchedScoringPlaneHits = new ArrayList<SimTrackerHit>();
        //Create a collection of LCRelations between track and scoring plane
        //hit
        List<LCRelation> trackToScoringplaneHitRelations = new ArrayList<LCRelation>();
        // Create a collection of LCRelations between a track and its corresponding MC particle
        List<LCRelation> trackToMCParticleRelations = new ArrayList<LCRelation>();

        MCParticle trackMCP = null;
        Map<Track, SimTrackerHit> trackScoringPlaneMap = new HashMap<Track,SimTrackerHit>(); 
        Map<Track, Cluster> trackClusterTruthMap = new HashMap<Track,Cluster>(); 
        Map<Track, MCParticle> trackMCParticleMap = new HashMap<Track, MCParticle>();
        Map<Cluster, MCParticle> clusterMCParticleMap = new HashMap<Cluster, MCParticle>();

        hitToRotated = TrackUtils.getHitToRotatedTable(event);
        hitToStrips = TrackUtils.getHitToStripsTable(event);
        List<TrackData> TrackData;
        RelationalTable TrktoData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
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

        List<MCParticle> genMCParticles = event.get(MCParticle.class, "MCParticle");
        for(MCParticle particle : genMCParticles){
            if(particle.getCharge() < 0)
                NgeneratedEle = NgeneratedEle + 1;
            else if (particle.getCharge() > 0)
                NgeneratedPos = NgeneratedPos + 1;
        }

        //Use matching algorithm to match tracks to clusters (without truth
        //info)

        Map<Track, Cluster> matchedTrackClusterMap = new HashMap<Track,Cluster>();
        matchedTrackClusterMap = matcher.newtrackClusterMatcher(tracks,TrktoData, hitToRotated, hitToStrips, this.trackCollectionName, clusters, this.trackClusterTimeOffset);

        for(Track track : tracks){
            
            int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
            if(charge < 0)
                NrecoEle = NrecoEle + 1;
            else
                NrecoPos = NrecoPos + 1;

            //Get track time
            double trackT;
            if (this.trackCollectionName.contains("GBLTracks")){
                trackT = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
            }
            else {
                trackdata = (TrackData) TrktoData.from(track);
                trackT = trackdata.getTrackTime();
            }



            /**
             * Get the MC particle associated with a track.
             * Fill mape with Track -> MCParticle
            **/
            trackMCP = this.getMCParticleAssociatedWithTrack(track,event);
            trackMCParticleMap.put(track, trackMCP);
            if(trackMCP == null) continue;


            // Add an LCRelation between the track and the corresponding MC particle
            trackToMCParticleRelations.add(new BaseLCRelation(track, trackMCP));

            /**
             * Match track to scoring plane hit if track_MCParticle ==
             * scoringPlaneHit_MCParticle.
             * scoring plane hit gives true track position at ECal.
             * */
            SimTrackerHit matchedScoringPlaneHit = null;
            for(SimTrackerHit scoringPlaneHit : scoringPlaneHits){
                // If the MC particles don't match, move on to the next particle
                if(!(scoringPlaneHit.getMCParticle() == trackMCP)) continue;
                // If a match is found, add the scoring plane hit to the list of matched hits and
                // an LCRelation between the track and the scoring plane.
                matchedScoringPlaneHits.add(scoringPlaneHit);
                trackToScoringplaneHitRelations.add(new BaseLCRelation(track, scoringPlaneHit));
                matchedScoringPlaneHit = scoringPlaneHit;
                //Build a map between Tracks and their scoringPlaneHits
                trackScoringPlaneMap.put(track, matchedScoringPlaneHit);

                // Once a match is found, there is no need to loop through the rest of the list
                break;
            }
            
            /**
             * Truth matching Tracks with Clusters via track_MCParticle ==
             * cluster_MCParticle.
             * */

            Cluster clusterTruthMatchedToTrack = null;
            for (Cluster cluster : clusters) {
                MCParticle clusterMCP = getMCParticlesAssociatedWithCluster(cluster, event);
                clusterMCParticleMap.put(cluster, clusterMCP);
                //if MCParticles of cluster and track match, map the two.
                if(clusterMCP == trackMCP) {
                    clusterTruthMatchedToTrack = cluster;                    
                    break;
                }
            }
            trackClusterTruthMap.put(track, clusterTruthMatchedToTrack);

            /**
             * Check purity and fake rate of KFTrackECalClusterMatcher
             * */
            if(clusterTruthMatchedToTrack != null){
                if(charge < 0)
                    NrecoElewMCP = NrecoElewMCP + 1;
                else
                    NrecoPoswMCP = NrecoPoswMCP + 1;

                trackClusterMatching(event, track, trackMCP, matchedTrackClusterMap.get(track), clusterTruthMatchedToTrack, clusterMCParticleMap.get(clusterTruthMatchedToTrack),clusters);
                trackClusterAnalysis(track,clusterTruthMatchedToTrack,trackT,"truth_matched");
            }
            else
                System.out.println("[EcalScoringPlaneDriver] NO TRUTH CLUSTER MATCHED TO TRACK!?");
        }
        trackScoringPlaneAnalysis(event,trackScoringPlaneMap);



        //event.put(ecalScoringPlaneHitsCollectionName, matchedScoringPlaneHits, SimTrackerHit.class, 0);
        //event.put(trackToScoringPlaneHitRelationsName, trackToScoringplaneHitRelations, LCRelation.class, 0);
        //event.put(trackToMCParticleRelationsName, trackToMCParticleRelations, LCRelation.class, 0);

    }

    public void trackClusterMatching(EventHeader event, Track track, MCParticle trackMCP,Cluster matchedCluster, Cluster truthmatchedCluster, MCParticle matchedclusterMCP, List<Cluster> clusters){
        
        //Get KFTrackData to access track time if KalmanFullTracks
        hitToRotated = TrackUtils.getHitToRotatedTable(event);
        hitToStrips = TrackUtils.getHitToStripsTable(event);
        List<TrackData> TrackData;
        RelationalTable TrktoData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
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
        //Get track parameters
        //double trackClusterTimeOffset = 44.8; //for 2016 MC
        int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
        double trackT;
        if (this.trackCollectionName.contains("GBLTracks")){
            trackT = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
        }
        else {
            trackdata = (TrackData) TrktoData.from(track);
            trackT = trackdata.getTrackTime();
        }

        //Map<Track,Cluster> trueMatches = new HashMap<Track,Cluster>();
        //Map<Track,Cluster> falseMatches = new HashMap<Track,Cluster>();
        //Cluster matchedCluster = matcher.trackClusterMatcher(track, this.trackCollectionName,charge,clusters, trackT, this.trackClusterTimeOffset);
        if(matchedCluster == null){
            System.out.println("[ScoringPlaneDriver] Failed to match recon track to EcalCluster");
            trackClusterAnalysis(track,matchedCluster,trackT,"no_match");
            if(charge > 0) 
                posNoMatch = posNoMatch + 1.0;
            else
                eleNoMatch = eleNoMatch + 1.0;
            return;
        }
        //MCParticle matchedClusterMCP = getMCParticlesAssociatedWithCluster(matchedCluster, event);;
                
        //if(trackMCP == matchedClusterMCP){
        if(matchedCluster == truthmatchedCluster){
            System.out.println("track_MCParticle matched to cluster_MCParticle");
            trackClusterAnalysis(track, matchedCluster,trackT,"positive_match");
            if(charge > 0)
                posTrueMatch = posTrueMatch + 1.0;
            else
                eleTrueMatch = eleTrueMatch + 1.0;
        }
        else{
            System.out.println("track_MCParticle failed to match to cluster_MCParticle");
            trackClusterAnalysis(track, matchedCluster,trackT,"negative_match");
            if(charge > 0) 
                posNegMatch = posNegMatch + 1.0;
            else
                eleNegMatch = eleNegMatch + 1.0;
        }

        eletotalCount = eleTrueMatch + eleNegMatch + eleNoMatch;
        postotalCount = posTrueMatch + posNegMatch + posNoMatch;
    }

    public void trackScoringPlaneAnalysis(EventHeader event, Map<Track,SimTrackerHit> trackScoringPlaneMap) {

        //Comparison plots for Track extrapolated to ECal vs Truth ScoringPlaneHit at Ecal
        hitToRotated = TrackUtils.getHitToRotatedTable(event);
        hitToStrips = TrackUtils.getHitToStripsTable(event);
        List<TrackData> TrackData;
        RelationalTable TrktoData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
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
    
        double trackx;
        double tracky;
        double trackz;
        double truthxpos;
        double truthypos;
        double truthzpos;
        double dxoffset;
        Map<Track, SimTrackerHit> map = trackScoringPlaneMap;

        for(Track track : map.keySet()) { 
            double trackT;
            double simTrackT;
            if (this.trackCollectionName.contains("GBLTracks")){
                trackT = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
            }
            else {
                trackdata = (TrackData) TrktoData.from(track);
                trackT = trackdata.getTrackTime();
            }
            SimTrackerHit matchedScoringPlaneHit = map.get(track);
            simTrackT = matchedScoringPlaneHit.getTime();

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

            truthxpos = matchedScoringPlaneHit.getPoint()[0];
            truthypos = matchedScoringPlaneHit.getPoint()[1];
            truthzpos = matchedScoringPlaneHit.getPoint()[2];

            double dx = truthxpos - trackx;
            double dy = truthypos - tracky;
            double dz = truthzpos - trackz;
            double dr = Math.sqrt(Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2));
            double dt = simTrackT - trackT;

            dxoffset = 0.0;

            //Make plots

            int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
            if(charge < 0) {


                //Track X,Y position plots at Ecal
                //
                plots2D.get(String.format("%s_ele_Track_Pos_at_Ecal; x mm; y mm",this.trackCollectionName)).fill(trackx,tracky);

                //Track at Ecal vs Cluster residuals
                plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dx",this.trackCollectionName)).fill(dx);
                plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dy",this.trackCollectionName)).fill(dy);
                plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dz",this.trackCollectionName)).fill(dz);
                plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dr",this.trackCollectionName)).fill(dr);
                plots1D.get(String.format("%s_ele_Track_atEcal_ScoringPlane_dt",this.trackCollectionName)).fill(dt);
            }
            else {
                //Track X,Y position at Ecal
                plots2D.get(String.format("%s_pos_Track_Pos_at_Ecal; x mm; y mm",this.trackCollectionName)).fill(trackx,tracky);

                //Track vs Cluster residuals at Ecal
                plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dx",this.trackCollectionName)).fill(dx);
                plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dy",this.trackCollectionName)).fill(dy);
                plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dz",this.trackCollectionName)).fill(dz);
                plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dr",this.trackCollectionName)).fill(dr);
                plots1D.get(String.format("%s_pos_Track_atEcal_ScoringPlane_dt",this.trackCollectionName)).fill(dt);
            }
        }

    }


    public void trackClusterAnalysis(Track track, Cluster cluster, double trackTime,  String identifier) {
        //For comparing extrapolated track position  with position of
        //truth-matched cluster
        String id = identifier; //positive_match, negative_match, truth_match
        double trackT = trackTime;
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
        //double[] trackP = new double[3];
        //double[] trackP = trackP = track.getMomentum();
        //System.out.println("trackP: " + trackP);
        
        if(trackCollectionName.contains("GBLTracks")) {
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
        //double trackPsum = Math.sqrt(Math.pow(trackP[0],2) + Math.pow(trackP[1],2) + Math.pow(trackP[2],2));

        dxoffset = 0.0;

        //if(cluster != null)
            //plots1D.get(String.format("%s_cluster_energy_%s",this.trackCollectionName,id)).fill(clusterEnergy);
        //plots1D.get(String.format("%s_track_momentum_%s",this.trackCollectionName,id)).fill(trackPsum);
        //plots1D.get(String.format("%s_momentum_energy_residual_%s",this.trackCollectionName,id)).fill(clusterEnergy - trackPsum);

        //Make plots
        int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
        if(id.contains("no_match")){
            /*
            if(charge < 0) {
                plots1D.get(String.format("%s_ele_track_NoMatchedCluster_x",this.trackCollectionName)).fill(trackx);
                plots1D.get(String.format("%s_ele_track_NoMatchedCluster_y",this.trackCollectionName)).fill(tracky);
                plots1D.get(String.format("%s_ele_track_NoMatchedCluster_z",this.trackCollectionName)).fill(trackz);
                plots1D.get(String.format("%s_ele_track_NoMatchedCluster_t",this.trackCollectionName)).fill(trackT);
            }
            else {

                plots1D.get(String.format("%s_pos_track_NoMatchedCluster_x",this.trackCollectionName)).fill(trackx);
                plots1D.get(String.format("%s_pos_track_NoMatchedCluster_y",this.trackCollectionName)).fill(tracky);
                plots1D.get(String.format("%s_pos_track_NoMatchedCluster_z",this.trackCollectionName)).fill(trackz);
                plots1D.get(String.format("%s_pos_track_NoMatchedCluster_t",this.trackCollectionName)).fill(dt);
            }
            */
            boolean skip = true;
        }
        else{
            if(charge < 0) {
                /*
                plots1D.get(String.format("%s_ele_track_MatchedCluster_x",this.trackCollectionName)).fill(trackx);
                plots1D.get(String.format("%s_ele_track_MatchedCluster_y",this.trackCollectionName)).fill(tracky);
                plots1D.get(String.format("%s_ele_track_MatchedCluster_z",this.trackCollectionName)).fill(trackz);
                plots1D.get(String.format("%s_ele_track_MatchedCluster_t",this.trackCollectionName)).fill(trackT);
                */
                plots1D.get(String.format("%s_ele_track_cluster_%s_dx",this.trackCollectionName,id)).fill(dx);
                plots1D.get(String.format("%s_ele_track_cluster_%s_dy",this.trackCollectionName,id)).fill(dy);
                plots1D.get(String.format("%s_ele_track_cluster_%s_dz",this.trackCollectionName,id)).fill(dz);
                plots1D.get(String.format("%s_ele_track_cluster_%s_dr",this.trackCollectionName,id)).fill(dr);
                plots1D.get(String.format("%s_ele_track_cluster_%s_dt",this.trackCollectionName,id)).fill(dt);
            }
            else {
                /*
                plots1D.get(String.format("%s_pos_track_MatchedCluster_x",this.trackCollectionName)).fill(trackx);
                plots1D.get(String.format("%s_pos_track_MatchedCluster_y",this.trackCollectionName)).fill(tracky);
                plots1D.get(String.format("%s_pos_track_MatchedCluster_z",this.trackCollectionName)).fill(trackz);
                plots1D.get(String.format("%s_pos_track_MatchedCluster_t",this.trackCollectionName)).fill(trackT);
                */
                plots1D.get(String.format("%s_pos_track_cluster_%s_dx",this.trackCollectionName,id)).fill(dx);
                plots1D.get(String.format("%s_pos_track_cluster_%s_dy",this.trackCollectionName,id)).fill(dy);
                plots1D.get(String.format("%s_pos_track_cluster_%s_dz",this.trackCollectionName,id)).fill(dz);
                plots1D.get(String.format("%s_pos_track_cluster_%s_dr",this.trackCollectionName,id)).fill(dr);
                plots1D.get(String.format("%s_pos_track_cluster_%s_dt",this.trackCollectionName,id)).fill(dt);
            }
        }


    }
/**
     * Get the MC particle associated with a track.
     * 
     * @param track : Track to get the MC particle for
     * @return The MC particle associated with the track
     */
    private MCParticle getMCParticlesAssociatedWithCluster(Cluster cluster, EventHeader event) {
        //Set<MCParticle>  ClusterUtilities.findMCParticles(cluster);
        //System.out.println("Finding MCParticles associated with Cluster");
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
        //List<MCParticle> mcParticles = new ArrayList<>();
        long cellID = seedhit.getCellID();
        for(RawTrackerHit ecalReadoutHit : ecalReadoutHits) {
            long recellID = ecalReadoutHit.getCellID();
            if(cellID == recellID) {
                readoutMatchHit = ecalReadoutHit;
                //System.out.println("CalorimeterHit " + cellID + " matched to ReadoutHit");
            }
        }


        Set<SimCalorimeterHit> simcalhits = rawtomc.allFrom(readoutMatchHit);
        double mcpEnergy = 0.0;
        MCParticle largestEnergyMCP = null;
        //System.out.println("Number of SimCalorimeterHits associated with readoutHit: " + simcalhits.size());
        for(SimCalorimeterHit simcalhit : simcalhits) {
            //System.out.println("MCParticle count per SimCalorimeterHit: " + simcalhit.getMCParticleCount());
            for(int i=0; i < simcalhit.getMCParticleCount(); i++){
                if(simcalhit.getMCParticle(i).getEnergy() > mcpEnergy) {
                    mcpEnergy = simcalhit.getMCParticle(i).getEnergy();
                    largestEnergyMCP = simcalhit.getMCParticle(i);
                    //mcParticles.add(simcalhit.getMCParticle(i));
                }
            }
        }
        //System.out.println("Cluster has " + mcParticles.size() + " MCParticles associated with it");

        return largestEnergyMCP;
    }
        

    private MCParticle getMCParticleAssociatedWithTrack(Track track, EventHeader event){
        
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
        System.out.println("[EcalScoringPlaneDriver]" + this.trackCollectionName + " NtrackerHits on Track: " + track.getTrackerHits().size());
        for(TrackerHit hit : track.getTrackerHits()){

            //int stripLayer = ((HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement()).getLayerNumber();
            int rawHitsLength = hit.getRawHits().size();
            List<RawTrackerHit> rawhits = hit.getRawHits();
            SimTrackerHit largestSimhit = null;
            double simhitMaxE = 0.0;
            for(RawTrackerHit rawhit : rawhits){
                Set<SimTrackerHit> simhits = rawtomc.allFrom(rawhit);
                for(SimTrackerHit simhit : simhits){
                    if (simhit != null && simhit.getMCParticle() != null) {
                        //System.out.println("simhit energy deposition: " + simhit.getdEdx());
                        double simhitEnergy = simhit.getdEdx();
                        if(simhitEnergy > simhitMaxE){
                            simhitMaxE = simhitEnergy;
                            largestSimhit = simhit;
                        }
                    }
                }
                            

                IDDecoder decoder = rawhit.getIDDecoder();
                decoder.setID(rawhit.getCellID());
                int layer = decoder.getLayer();
                SymmetricMatrix covmatrix = new SymmetricMatrix(3,hit.getCovMatrix(),false);
                BarrelEndcapFlag bef = BarrelEndcapFlag.UNKNOWN;
                BasicHep3Vector gpos = new BasicHep3Vector(hit.getPosition());
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
}


