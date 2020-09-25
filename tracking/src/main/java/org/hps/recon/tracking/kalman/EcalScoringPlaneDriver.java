//package org.hps.recon.test.ecalScoringPlane;
package org.hps.recon.tracking.kalman;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
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
import org.lcsim.event.GenericObject;


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
    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;
    KFTrackECalClusterMatcher matcher;

    boolean verbose = true;
    double elefakerate=0;
    double posfakerate=0;
    double elenegativeMatch=0;
    double elepositiveMatch=0;
    double posnegativeMatch=0;
    double pospositiveMatch=0;
    double eletotalCount=0;
    double postotalCount=0;

    //Collection Names
    String ecalScoringPlaneHitsCollectionName = "TrackerHitsECal";
    String tracksCollectionName = "KalmanFullTracks";
    String trackToScoringPlaneHitRelationsName = "TrackToEcalScoringPlaneHitRelations";
    String trackToMCParticleRelationsName = "TrackToMCParticleRelations";

    String ecalClustersCollectionName = "EcalClusters";
    String ecalReadoutHitsCollectionName = "EcalReadoutHits";
    String ecalTruthRelationsName = "EcalTruthRelations";

    private Set<SimTrackerHit> simhitsontrack = new HashSet<SimTrackerHit>();

    public void saveHistograms() {
        System.out.println("Saving Histogram");
        String rootFile = String.format("%s_TrackClusterMatching.root",this.tracksCollectionName);
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

        System.out.println("BOOKING HISTOGRAMS");
        String trackType = this.tracksCollectionName;
        plots1D = new HashMap<String, IHistogram1D>();
        tree = IAnalysisFactory.create().createTreeFactory().create();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

        plots1D.put(String.format("%s_ElectronTrackTruth_ECal_dx",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ElectronTrackTruth_ECal_dx",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_ElectronTrackTruth_ECal_dy",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ElectronTrackTruth_ECal_dy",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_ElectronTrackTruth_ECal_dz",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ElectronTrackTruth_ECal_dz",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_ElectronTrackTruth_ECal_dr",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ElectronTrackTruth_ECal_dr",this.tracksCollectionName), 100, -200, 200));

        plots1D.put(String.format("%s_PositronTrackTruth_ECal_dx",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_PositronTrackTruth_ECal_dx",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_PositronTrackTruth_ECal_dy",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_PositronTrackTruth_ECal_dy",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_PositronTrackTruth_ECal_dz",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_PositronTrackTruth_ECal_dz",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_PositronTrackTruth_ECal_dr",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_PositronTrackTruth_ECal_dr",this.tracksCollectionName), 100, -200, 200));


        plots1D.put(String.format("%s_ele_extrapTrack_truthCluster_dx",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_extrapTrack_truthCluster_dx",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_ele_extrapTrack_truthCluster_dy",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_extrapTrack_truthCluster_dy",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_ele_extrapTrack_truthCluster_dz",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_extrapTrack_truthCluster_dz",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_ele_extrapTrack_truthCluster_dr",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_extrapTrack_truthCluster_dr",this.tracksCollectionName), 100, -200, 200));

        plots1D.put(String.format("%s_pos_extrapTrack_truthCluster_dx",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_extrapTrack_truthCluster_dx",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_pos_extrapTrack_truthCluster_dy",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_extrapTrack_truthCluster_dy",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_pos_extrapTrack_truthCluster_dz",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_extrapTrack_truthCluster_dz",this.tracksCollectionName), 100, -200, 200));
        plots1D.put(String.format("%s_pos_extrapTrack_truthCluster_dr",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_extrapTrack_truthCluster_dr",this.tracksCollectionName), 100, -200, 200));


        plots1D.put(String.format("%s_pos_extrapTrack_truthCluster_dr",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_extrapTrack_truthCluster_dr",this.tracksCollectionName), 100, -200, 200));

        plots1D.put(String.format("%s_ele_fakeRate",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_fakeRate",this.tracksCollectionName), 10, -1, 10));
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void startOfData() {
        System.out.println("Starting job");
        bookHistograms();
        matcher = new KFTrackECalClusterMatcher(this.tracksCollectionName);
    }

    public void endOfData() {
        elefakerate = elenegativeMatch/eletotalCount;
        posfakerate = posnegativeMatch/postotalCount;

        for(int i=0; i < elefakerate; i++)
            plots1D.get(String.format("%s_ele_fakeRate",this.tracksCollectionName)).fill(1.0);
        for(int i=0; i < posfakerate; i++)
            plots1D.get(String.format("%s_ele_fakeRate",this.tracksCollectionName)).fill(1.0);
        saveHistograms();
    }

    protected void process(EventHeader event) {



        //Get collection of Ecal clusters
        //if (!event.hasCollection(Cluster.class, ecalClustersCollectionName)) return;
        //if (!event.hasCollection(Cluster.class, ecalTruthRelationsName)) return;
        //if (!event.hasCollection(Cluster.class, ecalReadoutHitsCollectionName)) return;

        List<Cluster> clusters = event.get(Cluster.class, ecalClustersCollectionName);

        //Find Track MCParticles and Scoring Plane hits
        
        //If event has no collection of tracks, skip
        if(!event.hasCollection(Track.class, tracksCollectionName)) return;
        //If even doesnt have collection of Ecal scoring plane hits, skip
        if(!event.hasCollection(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName)) return;

        // Get collection of tracks from event
        List<Track> tracks = event.get(Track.class, tracksCollectionName);
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

        for(Track track : tracks){

            // Get the MC particle associated with this track
            trackMCP = this.getMCParticleAssociatedWithTrack(track,event);
            trackMCParticleMap.put(track, trackMCP);
            SimTrackerHit matchedScoringPlaneHit = null;
            // If the MC particle is null, then the hits associated with the
            // track did not have an MC particle associated with them
            // TODO: Find out why some hits don't have any MC particles associated with them
            if(trackMCP == null) continue;

            // Add an LCRelation between the track and the corresponding MC particle
            trackToMCParticleRelations.add(new BaseLCRelation(track, trackMCP));

            // Loop over all of the scoring plane hits and check if the associated MC particle
            // matches the one from the track
            for(SimTrackerHit scoringPlaneHit : scoringPlaneHits){

                // If the MC particles don't match, move on to the next particle
                if(!(scoringPlaneHit.getMCParticle() == trackMCP)) continue;
                // If a match is found, add the scoring plane hit to the list of matched hits and
                // an LCRelation between the track and the scoring plane.
                matchedScoringPlaneHits.add(scoringPlaneHit);
                trackToScoringplaneHitRelations.add(new BaseLCRelation(track, scoringPlaneHit));
                matchedScoringPlaneHit = scoringPlaneHit;
                trackScoringPlaneMap.put(track, matchedScoringPlaneHit);

                // Once a match is found, there is no need to loop through the rest of the list
                break;
            }
            //Fill truth distance difference histogram
            //if(matchedScoringPlaneHit != null) {
              //  fillTrackScoringPlaneHisto(track, matchedScoringPlaneHit);
           // }
           
           

            //Truth Matching Tracks with Clusters

            List<Cluster> truthTrackMatchedCluster = new ArrayList<Cluster>();
            int index = 0;
            for (Cluster cluster : clusters) {
                MCParticle clusterMCP = getMCParticlesAssociatedWithCluster(cluster, event);
                clusterMCParticleMap.put(cluster, clusterMCP);
                if(clusterMCP == trackMCP) {
                    truthTrackMatchedCluster.add(cluster);                    
                    System.out.println("Matched truth track to truth cluster");
                    System.out.println("Cluster"+index+" MCParticle Charge: " + clusterMCP.getCharge());
                    System.out.println("Cluster"+index+" MCParticle Energy: " + clusterMCP.getEnergy());
                    //System.out.println("EndPosition: " + p.getEndPoint());
                    System.out.println("Cluster"+index+" Momentum: " + clusterMCP.getMomentum());

                    System.out.println("Track MCParticle Charge: " + trackMCP.getCharge());
                    System.out.println("Track MCParticle Energy: " + trackMCP.getEnergy());
                    //System.out.println("EndPosition: " + particle.getEndPoint());
                    System.out.println("Momentum: " + trackMCP.getMomentum());
                    break;
                }
                index = index + 1;
            }
            if(truthTrackMatchedCluster != null){
                index = 0;
                System.out.println("Length of truth matched Clusters: " + truthTrackMatchedCluster.size());
                for(Cluster cluster : truthTrackMatchedCluster) {
                    System.out.println("Cluster" + index + " Particle ID: " + cluster.getParticleId());
                    System.out.println("Cluster" + index + " Energy: " + cluster.getEnergy());
                    System.out.println("Cluster" + index + " Position: " + cluster.getPosition());
                    List<CalorimeterHit> calhits = cluster.getCalorimeterHits();
                    for(CalorimeterHit calhit : calhits)
                        System.out.println("Cluster" + index + "hit cellID's: " + calhit.getCellID());
                    index = index + 1;
                }
            }
            if(truthTrackMatchedCluster.size() > 0)
                trackClusterTruthMap.put(track, truthTrackMatchedCluster.get(0));
        }



        //event.put(ecalScoringPlaneHitsCollectionName, matchedScoringPlaneHits, SimTrackerHit.class, 0);
        //event.put(trackToScoringPlaneHitRelationsName, trackToScoringplaneHitRelations, LCRelation.class, 0);
        //event.put(trackToMCParticleRelationsName, trackToMCParticleRelations, LCRelation.class, 0);


        //ReconTrack to TruthTrack comparisons
        trackScoringPlaneAnalysis(trackScoringPlaneMap);
        //Compare extrapolated Track position at ECal to Truth-Matched ECal
        //cluster position
        reconTrackTruthClusterAnalysis(trackClusterTruthMap);
        //Check fake rate of track-cluster matching algorithm 
        trackClusterMatching(event, trackMCParticleMap, clusterMCParticleMap,clusters, tracks);

    }

    public void trackClusterMatching(EventHeader event, Map<Track, MCParticle> trackMCParticleMap, Map<Cluster, MCParticle> clusterMCParticleMap, List<Cluster> clusters, List<Track> tracks){
        
        hitToRotated = TrackUtils.getHitToRotatedTable(event);
        hitToStrips = TrackUtils.getHitToStripsTable(event);
        List<TrackData> tdata;
        RelationalTable TrktoData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
        List<LCRelation> trackRelations;
        GenericObject trackdata;
        tdata = event.get(TrackData.class, "KFTrackData");
        trackRelations = event.get(LCRelation.class, "KFTrackDataRelations");
        System.out.println("SIZE: " + trackRelations.size());
        for (LCRelation relation : trackRelations) {
            if (relation != null && relation.getTo() != null){
                TrktoData.add(relation.getFrom(), relation.getTo());
                System.out.println("LOOK: " + relation.getTo().toString());
                System.out.println("LOOK FROM: " + relation.getFrom().toString());
            }
        }


        double trackClusterTimeOffset = -43; //for 2016 MC
        for (Track track : tracks) {
            double trackT;
            int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());

            if (this.tracksCollectionName.contains("GBLTracks")){
                trackT = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
            }
            else {
                trackdata = (TrackData) TrktoData.from(track);
                //trackT = trackdata.getTrackTime();
                trackT = 0.0; 
            }

            Cluster matchedCluster = matcher.trackClusterMatcher(track, this.tracksCollectionName,charge,clusters, trackT, trackClusterTimeOffset);
            MCParticle matchedClusterMCP = null;
            for(Cluster cluster : clusterMCParticleMap.keySet()) {
                if(matchedCluster == cluster){
                    matchedClusterMCP = clusterMCParticleMap.get(cluster);
                    break;
                }
                    
            }
            if(trackMCParticleMap.get(track) == matchedClusterMCP){
                if(charge > 0)
                    pospositiveMatch = pospositiveMatch + 1.0;
                else
                    elepositiveMatch = elepositiveMatch + 1.0;
            }
            else{
                if(charge > 0) 
                    posnegativeMatch = posnegativeMatch + 1.0;
                else
                    elenegativeMatch = elenegativeMatch + 1.0;
            }

        }
        eletotalCount = elepositiveMatch + elenegativeMatch;
        postotalCount = pospositiveMatch + posnegativeMatch;

    }

    public void trackScoringPlaneAnalysis(Map<Track,SimTrackerHit> trackScoringPlaneMap) {

        //Comparison plots for Track extrapolated to ECal vs Truth ScoringPlaneHit at Ecal
    
        double trkxpos;
        double trkypos;
        double trkzpos;
        double truthxpos;
        double truthypos;
        double truthzpos;
        double dxoffset;
        Map<Track, SimTrackerHit> map = trackScoringPlaneMap;

        for(Track track : map.keySet()) { 
            SimTrackerHit matchedScoringPlaneHit = map.get(track);

            //Make histograms of truth vs extrapolation
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trkxpos = ts_ecalPos[0];
            trkypos = ts_ecalPos[1];
            trkzpos = ts_ecalPos[2];
            truthxpos = matchedScoringPlaneHit.getPoint()[0];
            truthypos = matchedScoringPlaneHit.getPoint()[1];
            truthzpos = matchedScoringPlaneHit.getPoint()[2];

            double dx = truthxpos - trkxpos;
            double dy = truthypos - trkypos;
            double dz = truthzpos - trkzpos;
            double dr = Math.sqrt(Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2));

            dxoffset = 0.0;

            //Make plots
            int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
            if(charge < 0) {
                plots1D.get(String.format("%s_ElectronTrackTruth_ECal_dx",this.tracksCollectionName)).fill(dx);
                plots1D.get(String.format("%s_ElectronTrackTruth_ECal_dy",this.tracksCollectionName)).fill(dy);
                plots1D.get(String.format("%s_ElectronTrackTruth_ECal_dz",this.tracksCollectionName)).fill(dz);
                plots1D.get(String.format("%s_ElectronTrackTruth_ECal_dr",this.tracksCollectionName)).fill(dr);
            }
            else {
                plots1D.get(String.format("%s_PositronTrackTruth_ECal_dx",this.tracksCollectionName)).fill(dx);
                plots1D.get(String.format("%s_PositronTrackTruth_ECal_dy",this.tracksCollectionName)).fill(dy);
                plots1D.get(String.format("%s_PositronTrackTruth_ECal_dz",this.tracksCollectionName)).fill(dz);
                plots1D.get(String.format("%s_PositronTrackTruth_ECal_dr",this.tracksCollectionName)).fill(dr);
            }
        }

    }


    public void reconTrackTruthClusterAnalysis(Map<Track, Cluster> trackClusterTruthMap) {
        //For comparing extrapolated track position  with position of
        //truth-matched cluster
        double trackx;
        double tracky;
        double trackz;
        double clustx;
        double clusty;
        double clustz;
        double dxoffset;
        Map<Track, Cluster> map = trackClusterTruthMap;
        
        for(Track track : map.keySet()) {
            Cluster cluster = map.get(track);
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trackx = ts_ecalPos[0];
            tracky = ts_ecalPos[1];
            trackz = ts_ecalPos[2];
            clustx = cluster.getPosition()[0];
            clusty = cluster.getPosition()[1];
            clustz = cluster.getPosition()[2];
            System.out.println("Cluster dx: " + clustx);
            System.out.println("Cluster dy: " + clusty);
            System.out.println("Cluster dz: " + clustz);

            double dx = trackx - clustx;
            double dy = tracky - clusty;
            double dz = trackz - clustz;
            double dr = Math.sqrt(Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2));
            System.out.println("LOOK HERE DX: " + dx);

            dxoffset = 0.0;

            //Make plots
            int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
            if(charge < 0) {
                plots1D.get(String.format("%s_ele_extrapTrack_truthCluster_dx",this.tracksCollectionName)).fill(dx);
                plots1D.get(String.format("%s_ele_extrapTrack_truthCluster_dy",this.tracksCollectionName)).fill(dy);
                plots1D.get(String.format("%s_ele_extrapTrack_truthCluster_dz",this.tracksCollectionName)).fill(dz);
                plots1D.get(String.format("%s_ele_extrapTrack_truthCluster_dr",this.tracksCollectionName)).fill(dr);
            }
            else {
                plots1D.get(String.format("%s_pos_extrapTrack_truthCluster_dx",this.tracksCollectionName)).fill(dx);
                plots1D.get(String.format("%s_pos_extrapTrack_truthCluster_dy",this.tracksCollectionName)).fill(dy);
                plots1D.get(String.format("%s_pos_extrapTrack_truthCluster_dz",this.tracksCollectionName)).fill(dz);
                plots1D.get(String.format("%s_pos_extrapTrack_truthCluster_dr",this.tracksCollectionName)).fill(dr);
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
                    System.out.println("rawtomc relations from: " + relation.getFrom().toString());
                    System.out.println("rawtomc relations to: " + relation.getTo().toString());
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
        for(TrackerHit hit : track.getTrackerHits()){

            //int stripLayer = ((HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement()).getLayerNumber();
            int rawHitsLength = hit.getRawHits().size();
            List<RawTrackerHit> rawhits = hit.getRawHits();
            for(RawTrackerHit rawhit : rawhits){
                Set<SimTrackerHit> simhits = rawtomc.allFrom(rawhit);
                IDDecoder decoder = rawhit.getIDDecoder();
                decoder.setID(rawhit.getCellID());
                int layer = decoder.getLayer();
                SymmetricMatrix covmatrix = new SymmetricMatrix(3,hit.getCovMatrix(),false);
                BarrelEndcapFlag bef = BarrelEndcapFlag.UNKNOWN;
                BasicHep3Vector gpos = new BasicHep3Vector(hit.getPosition());
                //HelicalTrackHit hth = new HelicalTrackHit(gpos, covmatrix, hit.getdEdx(), hit.getTime(), hit.getType(), hit.getRawHits(), "HPS-PhysicsRun2016-Pass2", layer, bef);
                for(SimTrackerHit simhit : simhits){
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
        return particle;
    }
}


