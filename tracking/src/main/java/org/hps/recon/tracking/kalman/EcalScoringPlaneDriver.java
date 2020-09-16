//package org.hps.recon.test.ecalScoringPlane;
package org.hps.recon.tracking.kalman;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;


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

        //Timing Plots
        plots1D.put(String.format("%s_ElectronTrack-Cluster_dt",this.tracksCollectionName), histogramFactory.createHistogram1D(String.format("%s_ElectronTrack-Cluster_dt",this.tracksCollectionName), 100, -200, 200));
    }
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void startOfData() {
        System.out.println("Starting job");
    }

    protected void process(EventHeader event) {
    
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
            double trkxpos;
            double trkypos;
            double trkzpos;
            double truthxpos;
            double truthypos;
            double truthzpos;
            double dxoffset;

            // Get the MC particle associated with this track
            particle = this.getMCParticleAssociatedWithTrack(track,event);
            SimTrackerHit matchedScoringPlaneHit = null;
            // If the MC particle is null, then the hits associated with the
            // track did not have an MC particle associated with them
            // TODO: Find out why some hits don't have any MC particles associated with them
            if(particle != null) System.out.println("Matching MC particle found!");
            if(particle == null) continue;

            // Add an LCRelation between the track and the corresponding MC particle
            trackToMCParticleRelations.add(new BaseLCRelation(track, particle));

            // Loop over all of the scoring plane hits and check if the associated MC particle
            // matches the one from the track
            for(SimTrackerHit scoringPlaneHit : scoringPlaneHits){

                // If the MC particles don't match, move on to the next particle
                if(!(scoringPlaneHit.getMCParticle() == particle)) continue;

                System.out.println("Found a match between a track and a scoring plane hit.");

                // If a match is found, add the scoring plane hit to the list of matched hits and
                // an LCRelation between the track and the scoring plane.
                matchedScoringPlaneHits.add(scoringPlaneHit);
                trackToScoringplaneHitRelations.add(new BaseLCRelation(track, scoringPlaneHit));
                matchedScoringPlaneHit = scoringPlaneHit;

                // Once a match is found, there is no need to loop through the rest of the list
                break;
            }

            //Make histograms of truth vs extrapolation
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trkxpos = ts_ecalPos[0];
            trkypos = ts_ecalPos[1];
            trkzpos = ts_ecalPos[2];
            truthxpos = matchedScoringPlaneHit.getPoint()[0];

            dxoffset = 0.0;

        }



        event.put(ecalScoringPlaneHitsCollectionName, matchedScoringPlaneHits, SimTrackerHit.class, 0);
        event.put(trackToScoringPlaneHitRelationsName, trackToScoringplaneHitRelations, LCRelation.class, 0);
        event.put(trackToMCParticleRelationsName, trackToMCParticleRelations, LCRelation.class, 0);
    }
/**
     * Get the MC particle associated with a track.
     * 
     * @param track : Track to get the MC particle for
     * @return The MC particle associated with the track
     */
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


