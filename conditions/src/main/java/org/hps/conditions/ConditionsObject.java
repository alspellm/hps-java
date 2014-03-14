package org.hps.conditions;

import java.sql.SQLException;

import org.hps.conditions.AbstractConditionsObject.FieldValueMap;

/**
 * This is an ORM interface for accessing conditions database information by row.
 * It can handle new or existing records.  The ID values for new records are -1
 * which indicates they are not in the database yet.
 * @author Jeremy McCormick <jeremym@slac.stanford.edu>
 */
public interface ConditionsObject {
    
    /**
     * Get the database table meta data associated to this object.
     * @return The database table meta data associated to this object.
     */
    ConditionsTableMetaData getTableMetaData();
    
    /**
     * Get the row ID of this object.
     * @return The database row ID.
     */
    int getRowId();
    
    /**
     * Get the collection ID of this object identifying its unique collection. 
     * @return The collection ID.
     */
    int getCollectionId();
    
    /**
     * Update this row in the database using a SQL UPDATE statement.    
     */
    void update() throws ConditionsObjectException;
    
    /**
     * Delete this object's row in the database using a SQL DELETE statement.
     */
    void delete() throws ConditionsObjectException;
    
    /**
     * Insert this object into the database using a SQL INSERT statement.     
     */
    void insert() throws ConditionsObjectException, SQLException;
    
    /**
     * Select data into this object from the database using a SQL SELECT statement.     
     */
    void select() throws ConditionsObjectException, SQLException;
    
    /**
     * Return true if this object is read-only.
     * @return True if object is read-only.
     */
    boolean isReadOnly();
    
    /**
     * Return true if this object is new and hasn't been inserted into the database yet.
     * @return True if object is new.
     */
    boolean isNew();
    
    /**
     * Return true if this object's data has been modified without a database update.
     * @return True if object is dirty.
     */
    boolean isDirty();
        
    /**
     * Generic set method for field values.  
     * This will set the object to the 'dirty' state.
     * @param fieldName The name of the field.
     * @param fieldValue The field value.
     */
    void setFieldValue(String field, Object value);
    
    /**
     * Set all of the field values on this object.
     * @param fieldValues The FieldValueMap containing pairs of names and values.
     */
    void setFieldValues(FieldValueMap fieldValues);
    
    /**
     * Get a field value, cast to the given type.
     * @param field The field value.
     * @return The field value casted to type T.
     */
    public <T> T getFieldValue(Class<T> klass, String field);
    
    /**
     * Set the ConnectionManager of this object.
     * This cannot be reset once set.
     * @param connectionManager The ConnectionManager.
     * @throws ConditionsObjectException if already set
     */
    void setConnectionManager(ConnectionManager connectionManager) throws ConditionsObjectException ;

    /**
     * Set the ConditionsTableMetaData of this object.
     * This cannot be reset once set.
     * @param tableMetaData The ConditionsTableMetaData.
     * @throws ConditionsObjectException if already set
     */
    void setTableMetaData(ConditionsTableMetaData tableMetaData) throws ConditionsObjectException;
    
    /**
     * Set the collection ID of this object.
     * This cannot be reset once set to a valid ID (e.g. not -1).
     * @param collectionId The collection ID.
     * @throws ConditionsObjectException if already set
     */
    void setCollectionId(int collectionId) throws ConditionsObjectException;
    
    /**
     * Set the row ID of this object.
     * This cannot be reset once set to a valid ID (e.g. not -1).
     * @param rowId The object's row ID.
     * @throws ConditionsObjectException if already set
     */
    public void setRowId(int rowId) throws ConditionsObjectException;
    
    /**
     * Set the object to read only mode.  This cannot be changed back once it is set.
     */
    void setIsReadOnly();
    
    /**
     * Generic Exception type throw by methods in this interface.
     */
    public static final class ConditionsObjectException extends Exception {
        public ConditionsObjectException(String message) {
            super(message);
        }
    }    
}
