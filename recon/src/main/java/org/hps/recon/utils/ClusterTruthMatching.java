package org.hps.recon.utils;

import java.util.ArrayList;
import java.util.HashMap;
//import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.lcsim.event.MCParticle;
import org.lcsim.event.RelationalTable;

import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.lcsim.event.Cluster;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.RawTrackerHit;

public class ClusterTruthMatching {

    boolean verbose = true;
    protected RawTrackerHit _readoutHit = null;
    protected SimCalorimeterHit _bestSimCalHit = null;
    protected MCParticle _bestMCP = null;
    private RelationalTable _rawtomc = null;

    protected Map<SimCalorimeterHit, List<MCParticle>> _simcalMCPMap = new HashMap<SimCalorimeterHit, List<MCParticle>>();
    protected Map<SimCalorimeterHit,Double> _seedhitSimhitsMap = new HashMap<SimCalorimeterHit, Double>();
    protected Map<MCParticle, Map<SimCalorimeterHit, Double>> _mcpEnergyContributionMap = new HashMap<MCParticle, Map<SimCalorimeterHit, Double>>();

    public ClusterTruthMatching(Cluster cluster, List<RawTrackerHit> ecalReadoutHits, RelationalTable rawtomc){
        this._rawtomc = rawtomc;
        doAnalysis(cluster, ecalReadoutHits, rawtomc);
    }

    private void doAnalysis(Cluster cluster, List<RawTrackerHit> ecalReadoutHits, RelationalTable rawtomc){
        getClusterReadoutHit(cluster, ecalReadoutHits); 
        getLargestSimCalHit(this._readoutHit, rawtomc);
        getLargestMcpFromSimCalHit(this._bestSimCalHit);
    }

    protected MCParticle getBestMCP(){
        return this._bestMCP;
    }

    protected MCParticle alternativeMethod(){
        getReadoutHitMCPs(this._readoutHit, this._rawtomc);
        MCParticle mcp = getLargestMCPEnergyContribution();
        return mcp;
    }   

    public void getClusterReadoutHit(Cluster cluster, List<RawTrackerHit> ecalReadoutHits){
        
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

        if(verbose){
            System.out.println("Checking Cluster with Energy: " + cluster.getEnergy());
            System.out.println("Cluster X: " + cluster.getPosition()[0]);
            System.out.println("Cluster Y: " + cluster.getPosition()[1]);
            System.out.println("Cluster seedhit cellID: " + cellID);
            System.out.println("Cluster seedhit Energy: " + seedhit.getRawEnergy());
        }

        if(readoutMatchHit == null){
            if(verbose)
                System.out.println("Error! No readouthit found for this cluster");
        }

        if(verbose){
            System.out.println("Matching Readout Hit Found: " + readoutMatchHit.getCellID());
            System.out.println("Readout Hit position: x= " + readoutMatchHit.getPosition()[0] + "; y= " + readoutMatchHit.getPosition()[1] + "; z= " + readoutMatchHit.getPosition()[2]);
        }
        this._readoutHit = readoutMatchHit;
    }

    public void getLargestSimCalHit(RawTrackerHit readoutHit, RelationalTable rawtomc){
        
        Set<SimCalorimeterHit> simcalhits = rawtomc.allFrom(readoutHit);
        double simcalhit_maxEnergy = 0.0;
        SimCalorimeterHit maxEnergy_simcalhit = null;
        double mcpMaxEnergy = 0.0;
        for(SimCalorimeterHit simcalhit : simcalhits){
            this._seedhitSimhitsMap.put(simcalhit, simcalhit.getRawEnergy());
            if(verbose)
                System.out.println("Simcalhit energy: " + simcalhit.getRawEnergy());
            if(simcalhit.getRawEnergy() > simcalhit_maxEnergy){
                simcalhit_maxEnergy = simcalhit.getRawEnergy();
                maxEnergy_simcalhit = simcalhit;
            }
        }

        if(verbose && maxEnergy_simcalhit != null){
            System.out.println("Simcalhit with energy: " + maxEnergy_simcalhit.getRawEnergy() + " selected");
            System.out.println("Simcalhit cellID: " + maxEnergy_simcalhit.getCellID());
            System.out.println("Simcalhit position: x= " + maxEnergy_simcalhit.getPosition()[0] + "; y= " + maxEnergy_simcalhit.getPosition()[1] + "; z= " + maxEnergy_simcalhit.getPosition()[2]);
        }

        this._bestSimCalHit = maxEnergy_simcalhit;
    }

    public void getReadoutHitMCPs(RawTrackerHit readoutHit, RelationalTable rawtomc){
        
        Set<SimCalorimeterHit> simcalhits = rawtomc.allFrom(readoutHit);
        for(SimCalorimeterHit simcalhit : simcalhits){
            List<MCParticle> mcps = new ArrayList<MCParticle>();
            for(int i=0; i < simcalhit.getMCParticleCount(); i++){
                MCParticle mcp = simcalhit.getMCParticle(i);
                mcps.add(mcp);
                double mcpEnergy = mcp.getEnergy();
                double mcpEnergyContr = simcalhit.getContributedEnergy(i);
                if(!_mcpEnergyContributionMap.containsKey(mcp)){
                    Map<SimCalorimeterHit, Double> tmp = new HashMap<SimCalorimeterHit, Double>();
                    tmp.put(simcalhit, mcpEnergyContr);
                    _mcpEnergyContributionMap.put(mcp, tmp);
                }
                else{
                    Map<SimCalorimeterHit, Double> tmp = _mcpEnergyContributionMap.get(mcp);
                    tmp.put(simcalhit, mcpEnergyContr);
                    _mcpEnergyContributionMap.put(mcp, tmp);
                }
            }       
            _simcalMCPMap.put(simcalhit, mcps);
        }
    }

    private MCParticle getLargestMCPEnergyContribution(){
        double maxEnergy = 0.0;
        MCParticle bestMCP = null;

        for(Map.Entry<MCParticle, Map<SimCalorimeterHit, Double>> entry : _mcpEnergyContributionMap.entrySet()){
            double mcpTotalE = 0.0;
            for(Map.Entry<SimCalorimeterHit, Double> subentry : entry.getValue().entrySet()){
                double energy = subentry.getValue();
                mcpTotalE = mcpTotalE + energy;
            }
            if(mcpTotalE > maxEnergy){
                maxEnergy = mcpTotalE;
                bestMCP = entry.getKey();
            }
        }
        return bestMCP;
    }

    public void getLargestMcpFromSimCalHit(SimCalorimeterHit simcalhit){
        
        List<MCParticle> mcps = new ArrayList<MCParticle>();

        double maxEnergy = 0.0;
        MCParticle bestMCP = null;
        for(int i=0; i < simcalhit.getMCParticleCount(); i++){
            MCParticle mcp = simcalhit.getMCParticle(i);
            double originZ = simcalhit.getMCParticle(i).getOriginZ();
            int pdgid = simcalhit.getMCParticle(i).getPDGID();
            double mcpEnergy = mcp.getEnergy();
            if(mcpEnergy > maxEnergy){
                maxEnergy = mcpEnergy;
                bestMCP = mcp;
            }
        }
        this._bestMCP = bestMCP;
    }
}








