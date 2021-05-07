package org.hps.readout.trigger2019;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import org.hps.readout.ReadoutDataManager;
import org.hps.readout.TriggerDriver;
import org.hps.readout.util.HodoscopePattern;
import org.hps.recon.ecal.EcalUtils;
import org.hps.record.daqconfig2019.ConfigurationManager2019;
import org.hps.record.daqconfig2019.DAQConfig2019;
import org.hps.record.triggerbank.TriggerModule2019;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.subdetector.HPSEcal3;

/**
 * <code>PairsTrigger2019ReadoutDriver</code> simulates an HPS pairs trigger
 * for 2019 MC. It takes in clusters produced by the
 * {@link org.hps.readout.ecal.updated.GTPClusterReadoutDriver
 * GTPClusterReadoutDriver}, and perform the necessary trigger
 * logic on them. If a trigger is detected, it is sent to the readout data
 * manager so that a triggered readout event may be written.
 */
public class PairsTrigger2019ReadoutDriver extends TriggerDriver{    
    
    /**
     * Indicates pair trigger type. Corresponding DAQ configuration is accessed by DAQ
     * configuration system, and applied into readout.
     */
    private String triggerType = "pair0";   
    
    // ==================================================================
    // ==== Trigger General Default Parameters ==========================
    // ==================================================================
    private String inputCollectionName = "EcalClustersGTP";       // Name for the LCIO cluster collection.
    private String inputCollectionNameHodo = "HodoscopePatterns";    
    
    private int pairCoincidence = 3;                              // Maximum allowed time difference between clusters. (4 ns clock-cycles)
    private String ecalGeometryName = "Ecal";                     // Name of the calorimeter geometry object.
    private TriggerModule2019 triggerModule = new TriggerModule2019();
    
    // ==================================================================
    // ==== Driver Internal Variables ===================================
    // ==================================================================
    private Queue<List<Cluster>> topClusterQueue = null;           // Store clusters on the top half of the calorimeter.
    private Queue<List<Cluster>> botClusterQueue = null;           // Store clusters on the bottom half of the calorimeter.
    private double localTime = 0.0;                                // Stores the internal time clock for the driver.
    private HPSEcal3 ecal = null;                                  // The calorimeter geometry object.
    
    @Override
    public void detectorChanged(Detector detector) {
        // Get the calorimeter sub-detector.
        org.lcsim.geometry.compact.Subdetector ecalSub = detector.getSubdetector(ecalGeometryName);
        if(ecalSub instanceof HPSEcal3) {
            ecal = (HPSEcal3) ecalSub;
        } else {
            throw new IllegalStateException("Error: Unexpected calorimeter sub-detector of type \"" + ecalSub.getClass().getSimpleName() + "; expected HPSEcal3.");
        }
    }
    
    /**
     * Sets whether or not the DAQ configuration is applied into the driver
     * the EvIO data stream or whether to read the configuration from data files.
     * 
     * @param state - <code>true</code> indicates that the DAQ configuration is
     * applied into the readout system, and <code>false</code> that it
     * is not applied into the readout system.
     */
    public void setDaqConfigurationAppliedintoReadout(boolean state) {
        // If the DAQ configuration should be read, attach a listener
        // to track when it updates.   
        if (state) {
            ConfigurationManager2019.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    // Get the DAQ configuration.
                    DAQConfig2019 daq = ConfigurationManager2019.getInstance();  
                    if(triggerType.contentEquals(PAIR0)) triggerModule.loadDAQConfiguration(daq.getVTPConfig().getPair0Config());
                    else if(triggerType.contentEquals(PAIR1)) triggerModule.loadDAQConfiguration(daq.getVTPConfig().getPair1Config());
                    else if(triggerType.contentEquals(PAIR2)) triggerModule.loadDAQConfiguration(daq.getVTPConfig().getPair2Config());
                    else if(triggerType.contentEquals(PAIR3)) triggerModule.loadDAQConfiguration(daq.getVTPConfig().getPair3Config());
                    pairCoincidence = (int)triggerModule.getCutValue(TriggerModule2019.PAIR_TIME_COINCIDENCE) / 4;
                }
            });
        }         
    }
    
    @Override
    public void startOfData() {                
        // Define the driver collection dependencies.
        addDependency(inputCollectionName);
        
        // Register the trigger.
        ReadoutDataManager.registerTrigger(this);
        
        // Make sure that a valid cluster collection name has been
        // defined. If it has not, throw an exception.
        if(inputCollectionName == null) {
            throw new RuntimeException("The parameter inputCollectionName was not set!");
        }
        
        // Initialize the top and bottom cluster queues.
        topClusterQueue = new LinkedList<List<Cluster>>();
        botClusterQueue = new LinkedList<List<Cluster>>();
        
        // Populate the top cluster queue. It should be populated with
        // a number of empty lists equal to (2 * pairCoincidence + 1).
        for(int i = 0; i < 2 * pairCoincidence + 1; i++) {
            topClusterQueue.add(new ArrayList<Cluster>());
        }
        
        // Populate the bottom cluster queue. It should be populated with
        // a number of empty lists equal to (2 * pairCoincidence + 1).
        for(int i = 0; i < pairCoincidence + 1; i++) {
            botClusterQueue.add(new ArrayList<Cluster>());
        }
        
        // Run the superclass method.
        super.startOfData();
    }
    
    @Override
    public void process(EventHeader event) {        
        // If there is no data ready, then nothing can be done/
        if(!ReadoutDataManager.checkCollectionStatus(inputCollectionName, localTime)) {
            return;
        }
        
        // Otherwise, get the input clusters from the present time.
        Collection<Cluster> clusters = ReadoutDataManager.getData(localTime, localTime + 4.0, inputCollectionName, Cluster.class);
        
        // Remove any clusters that do not pass the singles cuts. It
        // is more efficient to eliminate these before forming pairs.
        Collection<Cluster> goodClusters = getGoodClusters(clusters);
            
        // Put the good clusters into the cluster queue.
        updateClusterQueues(goodClusters);
        
        // Check that if a trigger exists, if the trigger is not in
        // dead time. If it is, no trigger may be issued, so this is
        // not necessary.
        if(!isInDeadTime() && testTrigger())
            sendTrigger(triggerType); 
        
        // Increment the local time.
        localTime += 4.0;
    }
    
    @Override
    protected double getTimeDisplacement() {
        return (pairCoincidence + 1) * 4.0;
    }
    
    @Override
    protected double getTimeNeededForLocalOutput() {
        return 0;
    }
    
    /**
     * Sets the name of the LCIO collection that contains the clusters.
     * @param clusterCollectionName - The cluster LCIO collection name.
     */
    public void setInputCollectionName(String clusterCollectionName) {
        inputCollectionName = clusterCollectionName;
    }
    
    public void setTriggerType(String trigger) {
        if(!trigger.equals(PAIR0) && !trigger.equals(PAIR1) && !trigger.equals(PAIR2) && !trigger.equals(PAIR3))
            throw new IllegalArgumentException("Error: wrong trigger type name \"" + trigger + "\".");
        triggerType = trigger;
    }
    
    /**
     * Sets the minimum number of hits needed for a cluster to pass
     * the hit count single cluster cut.
     * @param minHitCount - The parameter value.
     */
    public void setMinHitCount(int minHitCount) {
        triggerModule.setCutValue(TriggerModule2019.CLUSTER_HIT_COUNT_LOW, minHitCount);
    }
    
    /**
     * Sets the highest allowed energy a cluster may have and still
     * pass the cluster total energy single cluster cut. Value uses
     * units of GeV.
     * @param clusterEnergyHigh - The parameter value.
     */
    public void setClusterEnergyHigh(double clusterEnergyHigh) {
        triggerModule.setCutValue(TriggerModule2019.CLUSTER_TOTAL_ENERGY_HIGH, clusterEnergyHigh * EcalUtils.GeV);
    }
    
    /**
     * Sets the lowest allowed energy a cluster may have and still
     * pass the cluster total energy single cluster cut. Value uses
     * units of GeV.
     * @param clusterEnergyLow - The parameter value.
     */
    public void setClusterEnergyLow(double clusterEnergyLow) {
        triggerModule.setCutValue(TriggerModule2019.CLUSTER_TOTAL_ENERGY_LOW, clusterEnergyLow * EcalUtils.GeV);
    }
    
    /**
     * Sets the highest allowed energy a cluster pair may have and
     * still pass the cluster pair energy sum cluster cut. Value uses
     * units of GeV.
     * @param energySumHigh - The parameter value.
     */
    public void setEnergySumHigh(double energySumHigh) {
        triggerModule.setCutValue(TriggerModule2019.PAIR_ENERGY_SUM_HIGH, energySumHigh * EcalUtils.GeV);
    }
    
    /**
     * Sets the lowest allowed energy a cluster pair may have and
     * still pass the cluster pair energy sum cluster cut. Value uses
     * units of GeV.
     * @param energySumLow - The parameter value.
     */
    public void setEnergySumLow(double energySumLow) {
        triggerModule.setCutValue(TriggerModule2019.PAIR_ENERGY_SUM_LOW, energySumLow * EcalUtils.GeV);
    } 
    
    /**
     * Sets the highest allowed energy difference a cluster pair may
     * have and still pass the cluster pair energy difference cut.
     * Value uses units of GeV.
     * @param energyDifferenceHigh - The parameter value.
     */
    public void setEnergyDifferenceHigh(double energyDifferenceHigh) {
        triggerModule.setCutValue(TriggerModule2019.PAIR_ENERGY_DIFFERENCE_HIGH, energyDifferenceHigh * EcalUtils.GeV);
    }
    
    /**
     * Sets the lowest allowed energy slope a cluster pair may
     * have and still pass the cluster pair energy slope cut.
     * @param energySlopeLow - The parameter value.
     */
    public void setEnergySlopeLow(double energySlopeLow) {
        triggerModule.setCutValue(TriggerModule2019.PAIR_ENERGY_SLOPE_LOW, energySlopeLow);
    }
    
    /**
     * Sets the value of F in the energy slope equation <code>E_min +
     * R_min * F</code>.
     * @param f - The new energy slope parameter.
     */
    public void setEnergySlopeParamF(double f) {
        triggerModule.setCutValue(TriggerModule2019.PAIR_ENERGY_SLOPE_F, f);
    }
    
    /**
     * Sets the maximum deviation from coplanarity that a cluster pair
     * may possess and still pass the coplanarity pair cut. Value uses
     * units of degrees.
     * @param coplanarityHigh - The parameter value.
     */
    public void setCoplanarityHigh(double coplanarityHigh) {
        triggerModule.setCutValue(TriggerModule2019.PAIR_COPLANARITY_HIGH, coplanarityHigh);
    }

    /**
     * Sets the time range over which cluster pairs will be formed.
     * Value uses units of clock-cycles. Note that the number of
     * clock-cycles used is calculated as (2 * pairCoincidence) + 1.
     * @param pairCoincidence - The parameter value.
     */
    public void setPairCoincidence(int pairCoincidence) {
        this.pairCoincidence = pairCoincidence;
    }         
    
    private Collection<Cluster> getGoodClusters(Collection<Cluster> clusters) {
        // Create a list to hold clusters which pass the single
        // cluster cuts.
        Collection<Cluster> goodClusterList = new ArrayList<Cluster>();
        
        // Sort through the cluster list and add clusters that pass
        // the single cluster cuts to the good list.
        clusterLoop:
        for(Cluster cluster : clusters) {
            // ==== Cluster Hit Count Cut ==================================
            // =============================================================
            // If the cluster fails the cut, skip to the next cluster.
            if(triggerModule.getCutEn(TriggerModule2019.CLUSTER_HIT_COUNT_LOW_EN) && !triggerModule.clusterHitCountCut(cluster)) {
                continue clusterLoop;
            }
            
            // ==== Cluster Total Energy Cut ===============================
            // =============================================================
            // If the cluster fails the cut, skip to the next cluster.
            if(triggerModule.getCutEn(TriggerModule2019.CLUSTER_TOTAL_ENERGY_LOW_EN) && !triggerModule.clusterTotalEnergyCutLow(cluster)) {
                continue clusterLoop;
            }
            
            if(triggerModule.getCutEn(TriggerModule2019.CLUSTER_TOTAL_ENERGY_HIGH_EN) && !triggerModule.clusterTotalEnergyCutHigh(cluster)) {
                continue clusterLoop;
            }
            
            // A cluster that passes all of the single-cluster cuts
            // can be used in cluster pairs.
            goodClusterList.add(cluster);
        }
        
        // Return the good clusters.
        return goodClusterList;
    }
    
    /**
     * Adds clusters from a new event into the top and bottom cluster
     * queues so that they may be formed into pairs.
     * @param clusterList - The clusters to add to the queues.
     */
    private void updateClusterQueues(Collection<Cluster> clusterList) {
        // Create lists to store the top and bottom clusters.
        ArrayList<Cluster> topClusterList = new ArrayList<Cluster>();
        ArrayList<Cluster> botClusterList = new ArrayList<Cluster>();
        
        // Loop over the clusters in the cluster list.
        for(Cluster cluster : clusterList) {
            // Get the x and y indices. Note that LCSim meta data is
            // not available during readout, so crystal indices must
            // be obtained directly from the calorimeter geometry.
            java.awt.Point ixy = ecal.getCellIndices(cluster.getCalorimeterHits().get(0).getCellID());
            
            // If the cluster is on the top of the calorimeter, it
            // goes into the top cluster list.
            if(ixy.y > 0) { topClusterList.add(cluster); }
            
            // Otherwise, it goes into the bottom cluster list.
            else { botClusterList.add(cluster); }
        }
        
        // Add the new cluster lists to the cluster queues.
        topClusterQueue.add(topClusterList);
        botClusterQueue.add(botClusterList);
        
        // Remove the oldest cluster lists from the queues.
        topClusterQueue.remove();
        botClusterQueue.remove();
    }
    
    /**
     * Tests all of the current cluster pairs for triggers.
     * @return Returns <code>true</code> if one of the cluster pairs
     * passes all of the cluster cuts and <code>false</code> otherwise.
     */
    private boolean testTrigger() {
        // Track whether a trigger has occurred.
        boolean triggered = false;
        
        // Get the list of cluster pairs.
        List<Cluster[]> clusterPairs = getClusterPairsTopBot();
        
        // Iterate over the cluster pairs and perform each of the cluster
        // pair cuts on them. A cluster pair that passes all of the
        // cuts registers as a trigger.
        pairLoop:
        for(Cluster[] clusterPair : clusterPairs) {
            // Get the x and y indices. Note that LCSim meta data is
            // not available during readout, so crystal indices must
            // be obtained directly from the calorimeter geometry.
            java.awt.Point ixy0 = ecal.getCellIndices(clusterPair[0].getCalorimeterHits().get(0).getCellID());
            java.awt.Point ixy1 = ecal.getCellIndices(clusterPair[1].getCalorimeterHits().get(0).getCellID());
            
            
            // ==== Pair Energy Sum Cut ====================================
            // =============================================================
            // If the cluster fails the cut, skip to the next pair.
            if(triggerModule.getCutEn(TriggerModule2019.PAIR_ENERGY_SUM_EN) && !triggerModule.pairEnergySumCut(clusterPair)) {
                continue pairLoop;
            }
            
            // ==== Pair Energy Difference Cut =============================
            // =============================================================
            // If the cluster fails the cut, skip to the next pair.
            if(triggerModule.getCutEn(TriggerModule2019.PAIR_ENERGY_DIFFERENCE_HIGH_EN) && !triggerModule.pairEnergyDifferenceCut(clusterPair)) {
                continue pairLoop;
            }
            
            // ==== Pair Energy Slope Cut ==================================
            // =============================================================
            // If the cluster fails the cut, skip to the next pair.
            if(triggerModule.getCutEn(TriggerModule2019.PAIR_ENERGY_SLOPE_EN) && !triggerModule.pairEnergySlopeCut(clusterPair, ixy0, ixy1)) {
                continue pairLoop;
            }
            
            // ==== Pair Coplanarity Cut ===================================
            // =============================================================
            // If the cluster fails the cut, skip to the next pair.
            if(triggerModule.getCutEn(TriggerModule2019.PAIR_COPLANARITY_HIGH_EN) && !triggerModule.pairCoplanarityCut(clusterPair, ixy0, ixy1)) {
                continue pairLoop;
            }
            
            if(triggerType.equals(PAIR3)) {
                Collection<HodoscopePattern> hodoPatterns = null;
                ArrayList<HodoscopePattern> hodoPatternList = null;                
                if(ReadoutDataManager.checkCollectionStatus(inputCollectionNameHodo, localTime)) {
                    hodoPatterns = ReadoutDataManager.getData(localTime, localTime + 4.0, inputCollectionNameHodo, HodoscopePattern.class);
                    hodoPatternList = new ArrayList<>(hodoPatterns);  
                }
                                
                int clusterX0 = ixy0.x;
                if(clusterX0 < 0) clusterX0++;
                
                int clusterX1 = ixy1.x;
                if(clusterX1 < 0) clusterX1++;
                
                boolean cl0Stat = triggerModule.geometryMatchingCut(clusterX0, ixy0.y, hodoPatternList);
                boolean cl1Stat = triggerModule.geometryMatchingCut(clusterX1, ixy1.y, hodoPatternList);
                
                if(!cl0Stat && !cl1Stat) continue pairLoop;
            }
            
            // Clusters that pass all of the pair cuts produce a trigger.
            triggered = true;
        }
        
        // Return whether or not a trigger was observed.
        return triggered;
    }  
    
    /**
     * Get a list of all unique cluster pairs in the event.
     * @return A <code>List</code> collection of <code>Cluster</code> 
     * objects containing all cluster pairs.
     */
    protected List<Cluster[]> getClusterPairsTopBot() {
        // Create a list to store cluster pairs. 
        List<Cluster[]> clusterPairs = new ArrayList<Cluster[]>();
        
        // Loop over all top-bottom pairs of clusters; higher-energy cluster goes first in the pair
        // To apply pair coincidence time, use only bottom clusters from the 
        // readout cycle pairCoincidence readout cycles ago, and top clusters 
        // from all 2*pairCoincidence+1 previous readout cycles
        for(Cluster botCluster : botClusterQueue.element()) {
            for(List<Cluster> topClusters : topClusterQueue) {
                for(Cluster topCluster : topClusters) {
                    // The first cluster in a pair should always be
                    // the higher energy cluster. If the top cluster
                    // is higher energy, it goes first.
                    if(topCluster.getEnergy() > botCluster.getEnergy()) {
                        Cluster[] clusterPair = { topCluster, botCluster };
                        clusterPairs.add(clusterPair);
                    }
                    
                    // Otherwise, the bottom cluster goes first.
                    else {
                        Cluster[] clusterPair = { botCluster, topCluster };
                        clusterPairs.add(clusterPair);
                    }
                }
            }
        }
        
        // Return the cluster pair lists.
        return clusterPairs;
    }
    
}
