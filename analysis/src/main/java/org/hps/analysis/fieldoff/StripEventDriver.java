package org.hps.analysis.fieldoff;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TIData;
import org.lcsim.detector.DetectorElementStore;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.TrackerHit;
import org.lcsim.util.Driver;

/**
 *
 * @author Norman A Graf
 *
 */
public class StripEventDriver extends Driver {

    String _triggerType = "all";//allowed types are "" (blank) or "all", singles0, singles1, pairs0,pairs1

    private int _numberOfEventsWritten = 0;
    private int _minNumberOfStripHits = 0;
    private int _maxNumberOfStripHits = Integer.MAX_VALUE;
    private int _minNumberOfClusters = 0;
    private int _maxNumberOfClusters = Integer.MAX_VALUE;
    private double _minClusterEnergy = 0.;
    private double _maxClusterEnergy = 12.;
    private boolean _selectTopClusters = false;
    private boolean _selectBottomClusters = false;

    private String _clusterCollectionName = "EcalClustersCorr";

    private boolean _selectAllLayers = false;
    private int _requireNLayers = 12;
    private boolean _selectTopHits = false;
    private boolean _selectBottomHits = false;
    private boolean _requireTwoStripHitInLayerOne = false;

    @Override
    protected void process(EventHeader event) {
        boolean skipEvent = false;

        while (true) {
            if (!matchTrigger(event)) {
                skipEvent = true;
                break;
            }
            // let's start with clusters
            if (event.hasCollection(Cluster.class, _clusterCollectionName)) {
                List<Cluster> clusters = event.get(Cluster.class, _clusterCollectionName);
                int nclusters = clusters.size();
                if (nclusters < _minNumberOfClusters) {
                    skipEvent = true;
                    break;
                }
                if (nclusters > _maxNumberOfClusters) {
                    skipEvent = true;
                    break;
                }
                for (Cluster clus : clusters) {
                    double e = clus.getEnergy();
                    if (e < _minClusterEnergy) {
                        skipEvent = true;
                        break;
                    }
                    if (e > _maxClusterEnergy) {
                        skipEvent = true;
                        break;
                    }
                    double y = clus.getPosition()[1];
                    if (_selectTopClusters && y < 0) {
                        skipEvent = true;
                        break;
                    }
                    if (_selectBottomClusters && y > 0) {
                        skipEvent = true;
                        break;
                    }
                }
            }
            // now look at strip hits
            if (event.hasCollection(TrackerHit.class, "StripClusterer_SiTrackerHitStrip1D")) {
                int nHits = event.get(TrackerHit.class, "StripClusterer_SiTrackerHitStrip1D").size();
                if (nHits < _minNumberOfStripHits) {
                    skipEvent = true;
                    break;
                }
                if (nHits > _maxNumberOfStripHits) {
                    skipEvent = true;
                    break;
                }

                // add capability to require all tracking layers to have been hit
                if (_selectAllLayers) {
                    setupSensors(event);
                    Set<Integer> topLayers = new TreeSet<Integer>();
                    Set<Integer> bottomLayers = new TreeSet<Integer>();
                    int firstLayerAxialNstrips = 0;
                    int firstLayerStereoNstrips = 0;

                    List<TrackerHit> hits = event.get(TrackerHit.class, "StripClusterer_SiTrackerHitStrip1D");
                    for (TrackerHit hit : hits) {
                        List rawHits = hit.getRawHits();
                        HpsSiSensor sensor = null;
                        for (Object o : rawHits) {
                            RawTrackerHit rth = (RawTrackerHit) o;
                            sensor = (HpsSiSensor) rth.getDetectorElement();
                        }
                        String layerName = sensor.getName();
                        if (layerName.contains("t_")) {
                            topLayers.add(sensor.getLayerNumber());
                            // check on first layer
                            if (layerName.contains("L1t_")) {
                                if (layerName.contains("axial")) {
                                    firstLayerAxialNstrips = hit.getRawHits().size();
                                } else if (layerName.contains("stereo")) {
                                    firstLayerStereoNstrips = hit.getRawHits().size();
                                }
                            }
                        } else {
                            bottomLayers.add(sensor.getLayerNumber());
                            // check on first layer
                            if (layerName.contains("L1b_")) {
                                if (layerName.contains("axial")) {
                                    firstLayerAxialNstrips = hit.getRawHits().size();
                                } else if (layerName.contains("stereo")) {
                                    firstLayerStereoNstrips = hit.getRawHits().size();
                                }
                            }
                        }
                    }
                    if (_selectTopHits) {
                        if (topLayers.size() != _requireNLayers) {
                            skipEvent = true;
                            break;
                        }
                    }
                    if (_selectBottomHits) {
                        if (bottomLayers.size() != _requireNLayers) {
                            skipEvent = true;
                            break;
                        }
                    }
                    if (_requireTwoStripHitInLayerOne) {
                        if (firstLayerAxialNstrips != 2 || firstLayerStereoNstrips != 2) {
                            skipEvent = true;
                            break;
                        }
                    }
                    // if we don't explicitly request top or bottom, 
                    // only keep event if either the top or the bottom has all twelve layers hit
                    if (!_selectTopHits && !_selectBottomHits) {
                        skipEvent = true;
                        if (topLayers.size() == _requireNLayers) {
                            skipEvent = false;
                        }
                        if (bottomLayers.size() == _requireNLayers) {
                            skipEvent = false;
                        }
                    }
                }
            }
            // didn't fail any of the tests, so keep this event.
            skipEvent = false;
            break;
        }
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsWritten++;
        }
    }

    @Override
    protected void endOfData() {
        System.out.println("Wrote " + _numberOfEventsWritten + " events");
    }

    public void setMinNumberOfStripHits(int n) {
        _minNumberOfStripHits = n;
    }

    public void setMaxNumberOfStripHits(int n) {
        _maxNumberOfStripHits = n;
    }

    public void setMinNumberOfClusters(int n) {
        _minNumberOfClusters = n;
    }

    public void setMaxNumberOfClusters(int n) {
        _maxNumberOfClusters = n;
    }

    public void setMinClusterEnergy(double e) {
        _minClusterEnergy = e;
    }

    public void setMaxClusterEnergy(double e) {
        _maxClusterEnergy = e;
    }

    public void setClusterCollectionName(String s) {
        _clusterCollectionName = s;
    }

    public void setSelectTopHits(boolean b) {
        _selectTopHits = b;
    }

    public void setSelectBottomHits(boolean b) {
        _selectBottomHits = b;
    }

    public void setSelectTopClusters(boolean b) {
        _selectTopClusters = b;
    }

    public void setSelectBottomClusters(boolean b) {
        _selectBottomClusters = b;
    }

    public void setSelectAllLayers(boolean b) {
        _selectAllLayers = b;
    }

    public void setSelectNumberOfLayers(int i) {
        _requireNLayers = i;
    }

    public void setRequireTwoStripHitInLayerOne(boolean b) {
        _requireTwoStripHitInLayerOne = b;
    }

    private void setupSensors(EventHeader event) {
        List<RawTrackerHit> rawTrackerHits = event.get(RawTrackerHit.class, "SVTRawTrackerHits");
        EventHeader.LCMetaData meta = event.getMetaData(rawTrackerHits);
        // Get the ID dictionary and field information.
        IIdentifierDictionary dict = meta.getIDDecoder().getSubdetector().getDetectorElement().getIdentifierHelper().getIdentifierDictionary();
        int fieldIdx = dict.getFieldIndex("side");
        int sideIdx = dict.getFieldIndex("strip");
        for (RawTrackerHit hit : rawTrackerHits) {
            // The "side" and "strip" fields needs to be stripped from the ID for sensor lookup.
            IExpandedIdentifier expId = dict.unpack(hit.getIdentifier());
            expId.setValue(fieldIdx, 0);
            expId.setValue(sideIdx, 0);
            IIdentifier strippedId = dict.pack(expId);
            // Find the sensor DetectorElement.
            List<IDetectorElement> des = DetectorElementStore.getInstance().find(strippedId);
            if (des == null || des.size() == 0) {
                throw new RuntimeException("Failed to find any DetectorElements with stripped ID <0x" + Long.toHexString(strippedId.getValue()) + ">.");
            } else if (des.size() == 1) {
                hit.setDetectorElement((SiSensor) des.get(0));
            } else {
                // Use first sensor found, which should work unless there are sensors with duplicate IDs.
                for (IDetectorElement de : des) {
                    if (de instanceof SiSensor) {
                        hit.setDetectorElement((SiSensor) de);
                        break;
                    }
                }
            }
            // No sensor was found.
            if (hit.getDetectorElement() == null) {
                throw new RuntimeException("No sensor was found for hit with stripped ID <0x" + Long.toHexString(strippedId.getValue()) + ">.");
            }
        }
    }

    public void setTriggerType(String type) {
        _triggerType = type;
    }

    public boolean matchTriggerType(TIData triggerData) {
        if (_triggerType.contentEquals("") || _triggerType.contentEquals("all")) {
            return true;
        }
        if (triggerData.isSingle0Trigger() && _triggerType.contentEquals("singles0")) {
            return true;
        }
        if (triggerData.isSingle1Trigger() && _triggerType.contentEquals("singles1")) {
            return true;
        }
        if (triggerData.isPair0Trigger() && _triggerType.contentEquals("pairs0")) {
            return true;
        }
        if (triggerData.isPair1Trigger() && _triggerType.contentEquals("pairs1")) {
            return true;
        }
        if (triggerData.isPulserTrigger() && _triggerType.contentEquals("pulser")) {
            return true;
        }
        return false;

    }

    public boolean matchTrigger(EventHeader event) {
        boolean match = true;
        if (event.hasCollection(GenericObject.class, "TriggerBank")) {
            List<GenericObject> triggerList = event.get(GenericObject.class, "TriggerBank");
            for (GenericObject data : triggerList) {
                if (AbstractIntData.getTag(data) == TIData.BANK_TAG) {
                    TIData triggerData = new TIData(data);
                    if (!matchTriggerType(triggerData))//only process singles0 triggers...
                    {
                        match = false;
                    }
                }
            }
        }
        return match;
    }

}
