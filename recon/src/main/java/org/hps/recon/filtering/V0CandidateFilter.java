package org.hps.recon.filtering;

import static java.lang.Math.abs;
import java.util.List;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;

/**
 * Class to strip off trident candidates. Currently defined as: e+ e- events with
 * tracks matched to clusters. Neither electron can be a full-energy candidate
 * (momentum less than _fullEnergyCut [0.85GeV]). The Ecal cluster times must be 
 * within _timingCut [2.5ns] of each other.
 *
 * @author Norman A Graf
 *
 * @version $Id:
 */
public class V0CandidateFilter extends EventReconFilter {

    private String _V0CandidateCollectionName = "TargetConstrainedV0Candidates";
    private double _fullEnergyCut = 0.85;
    private double _clusterTimingCut = 2.5;

    private boolean _tight = false;

    @Override
    protected void process(EventHeader event) {
        incrementEventProcessed();
        if (!event.hasCollection(ReconstructedParticle.class, _V0CandidateCollectionName)) {
            skipEvent();
        }
        List<ReconstructedParticle> V0Candidates = event.get(ReconstructedParticle.class, _V0CandidateCollectionName);
        if (V0Candidates.size() == 0) {
            skipEvent();
        }

        // tight requires ONLY ONE real vertex fit 
        if (_tight) {
            if (V0Candidates.size() != 2) {
                skipEvent();
            }
        }

        for (ReconstructedParticle rp : V0Candidates) {

            ReconstructedParticle e1 = null;
            ReconstructedParticle e2 = null;

            List<ReconstructedParticle> electrons = rp.getParticles();
            if (electrons.size() != 2) {
                skipEvent();
            }
            // require both electrons to be associated with an ECal cluster
            e1 = electrons.get(0);
            if (e1.getClusters().size() == 0) {
                skipEvent();
            }
            e2 = electrons.get(1);
            if (e2.getClusters().size() == 0) {
                skipEvent();
            }
            // remove full energy electrons
            double p1 = e1.getMomentum().magnitude();
            if (p1 > _fullEnergyCut) {
                skipEvent();
            }
            double p2 = e2.getMomentum().magnitude();
            if (p2 > _fullEnergyCut) {
                skipEvent();
            }

            // calorimeter cluster timing cut
            // first CalorimeterHit in the list is the seed crystal
            double t1 = e1.getClusters().get(0).getCalorimeterHits().get(0).getTime();
            double t2 = e2.getClusters().get(0).getCalorimeterHits().get(0).getTime();

            if (abs(t1 - t2) > _clusterTimingCut) {
                skipEvent();
            }
        }
        incrementEventPassed();
    }

    /**
     * Maximum difference in Calorimeter Cluster Seed Hit times [ns]
     *
     * @param d
     */
    public void setClusterTimingCut(double d) {
        _clusterTimingCut = d;
    }

    /**
     * Name of V0 Candidate ReconstructedParticle Collection Name
     *
     * @param s
     */
    public void setV0CandidateCollectionName(String s) {
        _V0CandidateCollectionName = s;
    }

    /**
     * Maximum value for each of two electron momenta (removes full energy
     * electrons) [GeV]
     *
     * @param d
     */
    public void setV0MomentumMax(double d) {
        _fullEnergyCut = d;
    }

    /**
     * Setting a tight constraint requires one and only one candidate in the
     * event
     *
     * @param b
     */
    public void setTightConstraint(boolean b) {
        _tight = b;
    }
}