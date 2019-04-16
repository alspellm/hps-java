package org.hps.detector.hodoscope;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hps.conditions.hodoscope.HodoscopeChannel;
import org.hps.conditions.hodoscope.HodoscopeChannel.HodoscopeChannelCollection;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.IGeometryInfo;
import org.lcsim.detector.converter.compact.SubdetectorDetectorElement;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.identifier.IIdentifierHelper;
import org.lcsim.detector.solids.Box;

public final class HodoscopeDetectorElement extends SubdetectorDetectorElement{
                           
    private Map<Integer, List<HodoscopeChannel>> scintillatorPositionToChannelMap = new HashMap<Integer, List<HodoscopeChannel>>();    
    //private Map<Integer, String> hodoscopeBoundsMap = new HashMap<Integer, String>();
    
    public HodoscopeDetectorElement(String name, IDetectorElement parent) {
        super(name, parent);
    }
    
    /**
     * Gets the dimensions of the scintillator on which the hit
     * occurred. Returned values are one half the actual width,
     * height, and depth of the scintillator.
     * @param hit - The hit.
     * @return Returns a <code>double</code> array in the format of
     * <code>{ width / 2, height / 2, depth / 2 }</code>.
     */
    public final double[] getScintillatorHalfDimensions(IIdentifier id) {
        IDetectorElement idDetElem = findDetectorElement(id).get(0);
        Box box = (Box) idDetElem.getGeometry().getLogicalVolume().getSolid();
        return new double[] { box.getXHalfLength() + 0.05, box.getYHalfLength() + 0.05, box.getZHalfLength() + 0.25 };
    }
    
    /**
     * Gets the position of the center of the scintillator on which
     * the hit occurred.
     * @param hit - The hit.
     * @return Returns a <code>double</code> array in the format of
     * <code>{ x, y, z }</code>.
     */
    public final double[] getScintillatorPosition(IIdentifier id) {
        IDetectorElement idDetElem = findDetectorElement(id).get(0);
        IGeometryInfo geom = idDetElem.getGeometry();
        return new double[] { geom.getPosition().x(), geom.getPosition().y(), geom.getPosition().z() };
    }
    
    /**
     * Gets the x-, y-, and z-indices of the hodoscope scintillator
     * on which the hit occurred.
     * @param hit - The hit.
     * @return Returns the scintillator indices in the form of an
     * <code>int</code> array with the format <code>{ ix, iy, iz
     * }</code>.
     */
    public final int[] getHodoscopeIndices(IIdentifier id) {
        IIdentifierHelper helper = getIdentifierHelper();
        IIdentifierDictionary dict = helper.getIdentifierDictionary();
        IExpandedIdentifier expId = helper.unpack(id);
        return new int[] { expId.getValue(dict.getFieldIndex("ix")), 
                expId.getValue(dict.getFieldIndex("iy")), 
                expId.getValue(dict.getFieldIndex("layer")) };
    }
    
    /**
     * Gets the absolute positioning of the optical fiber hole(s) in
     * the scintillator in which a hit occurred. The array will be
     * either size 1 for scintillators with only one fiber hole or
     * size 2 for this with two. The first index corresponds to the
     * fiber hole that is closest to the center of the detector, and
     * the second to the one that is closer to the positron side of
     * the detector.
     * @param hit - The hit.
     * @return Returns a <code>double</code> array containing the
     * absolute positions of the optical fiber holes in the
     * scintillator.
     */
    public final double[] getScintillatorHolePositions(IIdentifier id) {
        // TODO: This is not ideal - ideally we would set this in the geometry somehow. This could hypothetically become invalid.
        // Get the x position at which the scintillator starts.
        double[] position = getScintillatorPosition(id);
        
        int ix = this.getIdentifierHelper().getValue(id, "ix");
        
        // The scintillators ix = 0 have only a single hole that is
        // in the exact geometric center of the scintillator.
        if(ix == 0) {
            return new double[] { position[0] };
        }
        
        // All other scintillators have two holes that all follow the
        // same paradigm.
        double[] halfDimension = getScintillatorHalfDimensions(id);
        double xStart = position[0] - halfDimension[0];
        double holeDisplacement = ((halfDimension[0] * 2) - 22.0) / 2;
        return new double[] { xStart + holeDisplacement, xStart + holeDisplacement + 22.0 };
    }
    
    /**
     * Creates a concise {@link java.lang.String String} object that
     * contains the geometric details for the scintillator on which
     * the argument hit occurred.
     * @param hit - The hit.
     * @return Returns a <code>String</code> describing the geometric
     * details of the scintillator on which the hit occurred.
     */
    public String getHodoscopeBounds(IIdentifier id) {
        int[] indices = getHodoscopeIndices(id);
        double[] dimensions = getScintillatorHalfDimensions(id);
        double[] position = getScintillatorPosition(id);
        double[] holeX = getScintillatorHolePositions(id);
        return String.format("Bounds for scintillator <%d, %2d, %d> :: x = [ %6.2f, %6.2f ], y = [ %6.2f, %6.2f ], z = [ %7.2f, %7.2f ];   FADC Channels :: %d;   Hole Positions:  x1 = %6.2f, x2 = %s",
                indices[0], indices[1], indices[2], position[0] - dimensions[0], position[0] + dimensions[0],
                position[1] - dimensions[1], position[1] + dimensions[1], position[2] - dimensions[2], position[2] + dimensions[2],
                getScintillatorChannelCount(id), holeX[0], holeX.length > 1 ? String.format("%6.2f", holeX[1]) : "N/A");
    }
    
    /**
     * Create a unique key for the hodoscope crystal.
     * @param ix - The x-index of the crystal. This should range from
     * 0 (for the closest to the beam) to 4 (for the farthest from
     * the beam).
     * @param iy - The y-index of the hodoscope. This may be either
     * -1 for the bottom half of the hodoscope or 1 for the top half.
     * @param layer - The layer number of the hodoscope. This may be
     * either 0 for layer 1 or 1 for layer 2.
     * @return Returns a unique integer representing these values.
     * The last bit represents the layer, the penultimate bit the top
     * or bottom position, and the remaining bits the x-index.
     */
    public static final int getHodoscopePositionVar(int ix, int iy, int layer) {
        int var = (ix << 2);
        var = var | ((iy == -1 ? 0 : 1) << 1);
        var = var | layer;
        return var;
    }
    
    /**
     * Gets the number of unique channels contained in the specified
     * scintillator. Some scintillators will have two optical fiber
     * holes which lead to separate channels. Others have either one
     * optical fiber hole or two, but each connect to the same FADC
     * channel.
     * @param ix - The x-index for the crystal.
     * @param iy - The y-index for the crystal.
     * @param iz - The layer number for the crystal.
     * @return Returns the number of unique FADC channels as an
     * <code>int</code>. This is <code>1</code> for scintillators
     * that have either one fiber hole or both fiber holes connect to
     * the same FADC channel. It is <code>2</code> otherwise.
     */
    public int getScintillatorChannelCount(int ix, int iy, int iz) {
        // Get the unique key for this scintillator.
        Integer posvar = Integer.valueOf(getHodoscopePositionVar(ix, iy, iz));
        
        // Return the number of unique channels.
        return scintillatorPositionToChannelMap.get(posvar).size();
    }
    
    /**
     * Gets the number of unique channels contained in the specified
     * scintillator. Some scintillators will have two optical fiber
     * holes which lead to separate channels. Others have either one
     * optical fiber hole or two, but each connect to the same FADC
     * channel.
     * @param hit - A hit occuring on the desired scintillator.
     * @return Returns the number of unique FADC channels as an
     * <code>int</code>. This is <code>1</code> for scintillators
     * that have either one fiber hole or both fiber holes connect to
     * the same FADC channel. It is <code>2</code> otherwise.
     */
    public int getScintillatorChannelCount(IIdentifier id) {
        int[] indices = getHodoscopeIndices(id);
        return getScintillatorChannelCount(indices[0], indices[1], indices[2]);
    }
            
    
    /**
     * Fills the mapping between scintillator position indices and
     * hardware FADC channel information.
     */
    public final void populateScintillatorChannelMap(HodoscopeChannelCollection channels) {
        for(HodoscopeChannel channel : channels) {
            Integer posvar = Integer.valueOf(getHodoscopePositionVar(channel.getX(), channel.getY(), channel.getLayer()));
            
            System.out.printf("Channel %2d at <%d, %2d, %d> and hole %d --> %s%n", channel.getChannelId(), channel.getX(), channel.getY(), channel.getLayer(), channel.getHole(), Integer.toBinaryString(posvar.intValue()));
            
            if(scintillatorPositionToChannelMap.containsKey(posvar)) {
                scintillatorPositionToChannelMap.get(posvar).add(channel);
            } else {
                List<HodoscopeChannel> channelList = new ArrayList<HodoscopeChannel>();
                channelList.add(channel);
                scintillatorPositionToChannelMap.put(posvar, channelList);
            }
        }
    }
    
    public Map<Integer, List<HodoscopeChannel>> getScintillatorPositionToChannelMap() {
        return this.scintillatorPositionToChannelMap;
    }
}