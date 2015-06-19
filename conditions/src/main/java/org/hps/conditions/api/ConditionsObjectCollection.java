package org.hps.conditions.api;

import java.util.Comparator;

/**
 * Interface representing a collection of conditions objects.
 *
 * @author Jeremy McCormick, SLAC
 * @param <ObjectType> the type of the objects
 */
public interface ConditionsObjectCollection<ObjectType extends ConditionsObject> extends Iterable<ObjectType>,
        DatabaseObject {

    /**
     * Add an object to the collection.
     *
     * @param object the object to add to the collection
     * @return <code>true</code> if object was added successfully
     * @throws ConditionsObjectException if there was an error adding the object
     */
    boolean add(final ObjectType object) throws ConditionsObjectException;

    /**
     * Add all objects to the collection.
     *
     * @param collection the source collection with objects to add
     */
    void addAll(ConditionsObjectCollection<ObjectType> collection);

    /**
     * Return <code>true</code> if collection contains this object.
     *
     * @param object the object to check
     * @return <code>true</code> if the collection contains the object
     */
    boolean contains(Object object);

    /**
     * Get an object by index.
     *
     * @param index the index of the object
     * @return the object
     */
    ObjectType get(final int index);

    /**
     * Get the collection ID.
     *
     * @return the collection ID
     */
    int getCollectionId();

    /**
     * Set the collection ID.
     *
     * @param collectionId the new collection ID
     */
    void setCollectionId(int collectionId);

    /**
     * Get the size of the collection.
     *
     * @return the size of the collection
     */
    int size();

    /**
     * Sort the collection in place.
     *
     * @param comparator the comparison operator to use for sorting
     */
    void sort(final Comparator<ObjectType> comparator);

    /**
     * Get a sorted copy of the collection, leaving the original in place.
     *
     * @param comparator the comparison operator to use
     * @return the sorted copy of the collection
     */
    ConditionsObjectCollection<ObjectType> sorted(final Comparator<ObjectType> comparator);
}