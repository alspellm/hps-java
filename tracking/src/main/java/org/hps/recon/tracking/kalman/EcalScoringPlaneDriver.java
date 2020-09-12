//package org.hps.recon.test.ecalScoringPlane;
package org.hps.recon.tracking.kalman;

//import java.util.HashMap;
import java.util.List; 
//import java.util.ArrayList; 
//import java.util.Map;

import org.lcsim.event.EventHeader;
//import org.lcsim.event.LCRelation;
//import org.lcsim.event.MCParticle;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Track;
//import org.lcsim.event.TrackerHit;
//import org.lcsim.event.base.BaseLCRelation;
//import org.lcsim.fit.helicaltrack.HelicalTrackHit;
import org.lcsim.util.Driver;

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
        
        if(verbose)  System.out.println("Track collection size: " + tracks.size());
        if(verbose)  System.out.println("Scoring Plane collection size: " + scoringPlaneHits.size());














    }
}
