//package org.hps.recon.test.ecalScoringPlane;
package org.hps.recon.tracking.kalman;

import java.util.HashMap;
import java.util.List; 
import java.util.ArrayList; 
import java.util.Map;

import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.base.BaseLCRelation;
import org.lcsim.fit.helicaltrack.HelicalTrackHit;
import org.lcsim.util.Driver;
//import org.lcsim.detector.IDetectorElement;
import org.lcsim.geometry.IDDecoder;

//import hep.physics.vec.BasicHep3Vector;
//import hep.physics.vec.Hep3Vector;

//import org.hps.recon.tracking.TrackUtils;

/** 
 * Driver stolen from Omar to relate a Track to an Ecal scoring plane hit
 *
 **/

public class EcalScoringPlaneDriver extends Driver {

    boolean verbose = true;

    //Collection Names
    String ecalScoringPlaneHitsCollectionName = "TrackerHitsECal";
    String tracksCollectionName = "MatchedTracks";
    String trackToScoringPlaneHitRelationsName = "TrackToEcalScoringPlaneHitRelations";
    String trackToMCParticleRelationsName = "TrackToMCParticleRelations";

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

            // Get the MC particle associated with this track
            particle = this.getMCParticleAssociatedWithTrack(track);
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

                // Once a match is found, there is no need to loop through the rest of the list
                break;
            }
        }


        
        if(verbose)  System.out.println("Track collection size: " + tracks.size());
        if(verbose)  System.out.println("Scoring Plane collection size: " + scoringPlaneHits.size());

    }
/**
     * Get the MC particle associated with a track.
     * 
     * @param track : Track to get the MC particle for
     * @return The MC particle associated with the track
     */
    private MCParticle getMCParticleAssociatedWithTrack(Track track){
        
        Map <MCParticle, int[]>mcParticleMultiplicity = new HashMap<MCParticle, int[]>();
        MCParticle particle;
        for(TrackerHit hit : track.getTrackerHits()){


            //get raw hit -> get Sim tracker hit?
            //HelicalTrackHit hit = TrackUtils.makeHelicalTrackHitFromTrackerHit(trhit);
            // If one of the tracker hits doesn't have any MC particles associated
            // with it, return null for now.
            //HelicalTrackHit hth = new HelicalTrackHit(hit.getPosition(),
            //hit.getCovMatrix(), hit.getdEdx, hit.getTime(), hit.getType(), hit.getRawHits(), "HPS-PhysicsRun2016-Pass2");
            // BasicHep3Vector gpos = new BasicHep3Vector(hit.getPosition());
            List<RawTrackerHit> rawhits = hit.getRawHits();
            //System.out.println("TYPE " + rawhits.get(0).getType());
            //System.out.println("DE " + hit.getRawHits().get(0).getDetectorElement());

            //System.out.println("Detector Element: " + IDetectorElement.findDetectorElement((Hep3Vector) gpos));
            if(((HelicalTrackHit) hit).getMCParticles().size() == 0){
                System.out.println("HelicalTrackHit is not associated with any MC particles.");
                return null;
            }
            else System.out.println("GOT ONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            
            particle = ((HelicalTrackHit) hit).getMCParticles().get(0);
            if(!mcParticleMultiplicity.containsKey(particle)){
                mcParticleMultiplicity.put(particle, new int[1]);
                mcParticleMultiplicity.get(particle)[0] = 0;
            }
            
            mcParticleMultiplicity.get(particle)[0]++;
                
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


