//package org.hps.recon.test.ecalScoringPlane;
package org.hps.recon.tracking.kalman;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;

//import org.hps.recon.ecal.cluster.ClusterUtilities;


import java.util.HashMap;
import java.util.List; 
import java.util.ArrayList; 
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.io.IOException;

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

    boolean verbose = true;

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
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void startOfData() {
        System.out.println("Starting job");
        bookHistograms();
    }

    public void endOfData() {
        saveHistograms();
    }

    protected void process(EventHeader event) {
    

        //Get collection of Ecal clusters
        //if (!event.hasCollection(Cluster.class, ecalClustersCollectionName)) return;
        //if (!event.hasCollection(Cluster.class, ecalTruthRelationsName)) return;
        //if (!event.hasCollection(Cluster.class, ecalReadoutHitsCollectionName)) return;

        List<Cluster> clusters = event.get(Cluster.class, ecalClustersCollectionName);
        Map<Cluster, List<MCParticle>> clusterMCMap = new HashMap<Cluster,List<MCParticle>>(); 

        System.out.println("Number of clusters in event:" + clusters.size());
        /*
        for (Cluster cluster : clusters) {
            List<MCParticle> mcParticles = getMCParticlesAssociatedWithCluster(cluster, event);
            System.out.println("Number of MCParticles associated with cluster: " + mcParticles.size()); 
            for(MCParticle mcParticle : mcParticles) {
                if (clusterMCMap.get(cluster) == null) {
                    clusterMCMap.put(cluster, new ArrayList<MCParticle>());
                }
                clusterMCMap.get(cluster).add(mcParticle);

            }
        }
        */
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
        MCParticle particle = null;
        for(Track track : tracks){

            // Get the MC particle associated with this track
            particle = this.getMCParticleAssociatedWithTrack(track,event);
            SimTrackerHit matchedScoringPlaneHit = null;
            // If the MC particle is null, then the hits associated with the
            // track did not have an MC particle associated with them
            // TODO: Find out why some hits don't have any MC particles associated with them
            if(particle == null) continue;

            // Add an LCRelation between the track and the corresponding MC particle
            trackToMCParticleRelations.add(new BaseLCRelation(track, particle));

            // Loop over all of the scoring plane hits and check if the associated MC particle
            // matches the one from the track
            for(SimTrackerHit scoringPlaneHit : scoringPlaneHits){

                // If the MC particles don't match, move on to the next particle
                if(!(scoringPlaneHit.getMCParticle() == particle)) continue;
                // If a match is found, add the scoring plane hit to the list of matched hits and
                // an LCRelation between the track and the scoring plane.
                matchedScoringPlaneHits.add(scoringPlaneHit);
                trackToScoringplaneHitRelations.add(new BaseLCRelation(track, scoringPlaneHit));
                matchedScoringPlaneHit = scoringPlaneHit;

                // Once a match is found, there is no need to loop through the rest of the list
                break;
            }
            //Fill truth distance difference histogram
            if(matchedScoringPlaneHit != null) {
                fillTruthDistanceHisto(track, matchedScoringPlaneHit);
            }

            List<Cluster> truthTrackClusterMatch = new ArrayList<Cluster>();
            int index = 0;
            for (Cluster cluster : clusters) {
                List<MCParticle> mcParticles = getMCParticlesAssociatedWithCluster(cluster, event);
                for(MCParticle mcparticle : mcParticles) {
                    if(mcparticle == particle) {
                        truthTrackClusterMatch.add(cluster);                    
                        System.out.println("Matched truth track to truth cluster");
                        System.out.println("Cluster"+index+" MCParticle Charge: " + mcparticle.getCharge());
                        System.out.println("Cluster"+index+" MCParticle Energy: " + mcparticle.getEnergy());
                        //System.out.println("EndPosition: " + p.getEndPoint());
                        System.out.println("Cluster"+index+" Momentum: " + mcparticle.getMomentum());

                        System.out.println("Track MCParticle Charge: " + particle.getCharge());
                        System.out.println("Track MCParticle Energy: " + particle.getEnergy());
                        //System.out.println("EndPosition: " + particle.getEndPoint());
                        System.out.println("Momentum: " + particle.getMomentum());
                        break;
                    }
                }
                index = index + 1;

            }
            if(truthTrackClusterMatch != null){
                index = 0;
                System.out.println("Length of truth matched Clusters: " + truthTrackClusterMatch.size());
                for(Cluster cluster : truthTrackClusterMatch) {
                    System.out.println("Cluster" + index + " Particle ID: " + cluster.getParticleId());
                    System.out.println("Cluster" + index + " Energy: " + cluster.getEnergy());
                    System.out.println("Cluster" + index + " Position: " + cluster.getPosition());
                    List<CalorimeterHit> calhits = cluster.getCalorimeterHits();
                    for(CalorimeterHit calhit : calhits)
                        System.out.println("Cluster" + index + "hit cellID's: " + calhit.getCellID());
                    index = index + 1;
                }
            }


            /*    
            //Use Map of Cluster -> MCParticles and MCParticle -> Track to
            //match Track and Cluster
            int repeatCluster=-1; 
            for(Map.Entry<Cluster, List<MCParticle>> entry : clusterMCMap.entrySet()) {
                for(MCParticle p : entry.getValue()){
                    System.out.println("Clusters inside of Cluster: " + entry.getKey().getClusters().size());
                    if(p.getCharge() == particle.getCharge()) {
                        if(p == particle) {
                            System.out.println("Matching Cluster to Track MCParticle Found");
                            System.out.println("Cluster MCParticle Charge: " + p.getCharge());
                            System.out.println("Cluster MCParticle Energy: " + p.getEnergy());
                            //System.out.println("EndPosition: " + p.getEndPoint());
                            System.out.println("Cluster Momentum: " + p.getMomentum());
                        }
                    }
                    
                }
            }
            /*
            for(Map.Entry<Cluster, List<MCParticle>> entry : clusterMCMap.entrySet()) {
                for(MCParticle p : entry.getValue())
                {
                    if(p == particle) {
                        System.out.println("Matched Track MCParticle to Cluster MCParticle");
                        repeatCluster = repeatCluster + 1;
                        break;
                    }
                }
            }
            if(repeatCluster > 0)
                System.out.println("Multiple Clusters Matched to Track. No Bueno");
                */
        }



        event.put(ecalScoringPlaneHitsCollectionName, matchedScoringPlaneHits, SimTrackerHit.class, 0);
        event.put(trackToScoringPlaneHitRelationsName, trackToScoringplaneHitRelations, LCRelation.class, 0);
        event.put(trackToMCParticleRelationsName, trackToMCParticleRelations, LCRelation.class, 0);


    }

    public void fillTruthDistanceHisto(Track trk, SimTrackerHit scoringPlaneHit) {
    
        double trkxpos;
        double trkypos;
        double trkzpos;
        double truthxpos;
        double truthypos;
        double truthzpos;
        double dxoffset;
        Track track = trk;
        SimTrackerHit matchedScoringPlaneHit = scoringPlaneHit;

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
/**
     * Get the MC particle associated with a track.
     * 
     * @param track : Track to get the MC particle for
     * @return The MC particle associated with the track
     */
    private List<MCParticle> getMCParticlesAssociatedWithCluster(Cluster cluster, EventHeader event) {
        //Set<MCParticle>  ClusterUtilities.findMCParticles(cluster);
        //System.out.println("Finding MCParticles associated with Cluster");
        List <RawTrackerHit> ecalReadoutHits = event.get(RawTrackerHit.class,ecalReadoutHitsCollectionName);
        Map <MCParticle, int[]>mcParticleMultiplicity = new HashMap<MCParticle, int[]>();
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, ecalTruthRelationsName )) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, ecalTruthRelationsName);
            for (LCRelation relation : trueHitRelations)
                if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                    rawtomc.add(relation.getFrom(), relation.getTo());
                else System.out.println("failed to build cluster mc relation");
        }

        //Match Cluster->CalorimeterHit(s) to corresponding EcalReadoutHits via
        //CellID match

        List<CalorimeterHit> calorimeterhits = cluster.getCalorimeterHits();
        List<MCParticle> mcParticles = new ArrayList<>();
        for(CalorimeterHit calorimeterhit : calorimeterhits) {
            RawTrackerHit readoutMatchHit = null;
            int repeat=-1;
            long cellID = calorimeterhit.getCellID();
            for(RawTrackerHit ecalReadoutHit : ecalReadoutHits) {
                long recellID = ecalReadoutHit.getCellID();
                if(cellID == recellID) {
                    readoutMatchHit = ecalReadoutHit;
                    repeat = repeat + 1;
                    //System.out.println("CalorimeterHit " + cellID + " matched to ReadoutHit");
                }
            }
            if(repeat > 0){
                System.out.println("Bad Cluster Readout Match");
                return null;
            }

            Set<SimCalorimeterHit> simcalhits = rawtomc.allFrom(readoutMatchHit);
            //System.out.println("Number of SimCalorimeterHits associated with readoutHit: " + simcalhits.size());
            for(SimCalorimeterHit simcalhit : simcalhits) {
                //System.out.println("MCParticle count per SimCalorimeterHit: " + simcalhit.getMCParticleCount());
                for(int i=0; i < simcalhit.getMCParticleCount(); i++){
                    mcParticles.add(simcalhit.getMCParticle(i));
                    //System.out.println("MCP Momentum: " + simcalhit.getMCParticle(i).getMomentumAtEndpoint());
                    //System.out.println("MCP Endpoint: " + simcalhit.getMCParticle(i).getEndPoint());
                    //System.out.println("MCP Charge: " + simcalhit.getMCParticle(i).getCharge());
                    //System.out.println("MCP PDGID: " + simcalhit.getMCParticle(i).getPDGID());
                }
            }
        }
        //System.out.println("Cluster has " + mcParticles.size() + " MCParticles associated with it");

        return mcParticles;
        

        /*
        for(RawTrackerHit rawCalHit : rawCalHits) {
            System.out.println("got a rawCalHit");
             
        }

        CalorimeterHit seedhit = ClusterUtilities.findSeedHit(cluster);

        System.out.println("Running getMCParticleAssociatedWithCluster");
        for(CalorimeterHit calorimeterhit : calorimeterhits){
            double hitrawenergy = calorimeterhit.getCorrectedEnergy();
            System.out.println("hit energy: " + hitrawenergy);
            for(SimCalorimeterHit simCalHit : simCalHits){
                double simrawenergy = simCalHit.getCorrectedEnergy();
                System.out.println("sim energy: " + simrawenergy);
                if(hitrawenergy == simrawenergy) System.out.println("Matching Energies Found!");
            }
        }

        System.out.println("seedhit id: " + seedhit.getCellID());
        RawCalorimeterHit rawseedhit = (RawCalorimeterHit) seedhit;
        if(simcalhits != null) {
        Set<SimCalorimeterHit> simcalhits = rawtomc.allFrom(rawseedhit);
            System.out.println("simcalhits found"); 
        }
        for(SimCalorimeterHit simcalhit : simcalhits) {
            System.out.println("Cluster MCParticle count: " + simcalhit.getMCParticleCount());
            MCParticle mcparticle = simcalhit.getMCParticle(1);
            if(mcparticle != null) System.out.println("Found cluster MCParticle!");
        }
        Set<MCParticle> particles = ClusterUtilities.findMCParticles(cluster);
        for(MCParticle particle : particles) {
            if(particle != null) System.out.println("FINALLY FOUND AN MCPARTICLE");
        }

        for(CalorimeterHit calorimeterhit : calorimeterhits) {
            RawCalorimeterHit rawcalorimeterhit = (RawCalorimeterHit) calorimeterhit;
            Set<SimCalorimeterHit> simcalhits = rawtomc.allFrom(rawcalorimeterhit);
            System.out.println("simcalhits size: " + simcalhits.size());
            for (SimCalorimeterHit simcalhit : simcalhits) {
                System.out.println("Cluster simcalhit found");
                if (simcalhit != null && simcalhit.getMCParticle(1) != null){
                    System.out.println("Cluster MCParticle found!");
                    break;
                }
            }
        }
            */

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
            System.out.println("size of raw hits: " + rawHitsLength);
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
                        System.out.println("FOUND A SIMHIT!");
                        particle = simhit.getMCParticle();
                        if(particle != null) System.out.println("FOUND A MATCHING MCPARTICLE");
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


