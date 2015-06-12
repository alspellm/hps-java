package org.hps.conditions.api;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.SimpleDateFormat;
import java.util.Date;

import org.hps.conditions.database.Field;

/**
 * This is a basic ORM class for performing CRUD (create, read, update, delete) operations on objects in the conditions
 * system. Each object is mapped to a single row in a database table.
 *
 * @author <a href="mailto:jeremym@slac.stanford.edu">Jeremy McCormick</a>
 */
public class BaseConditionsObject implements ConditionsObject {

    /**
     * Field name for collection ID.
     */
    static final String COLLECTION_ID_FIELD = "collection_id";

    /**
     * Date formatting for insert statement.
     */
    static final SimpleDateFormat DATE_FORMAT = new SimpleDateFormat("yyyy-MM-dd kk:mm:ss");

    /**
     * Value that indicates collection ID is not set (no collection assigned).
     */
    static final int UNSET_COLLECTION_ID = -1;

    /**
     * Value that indicates row ID is not assigned (new record).
     */
    static final int UNSET_ROW_ID = -1;

    /**
     * Perform the default to string operation on a conditions object.
     * 
     * @param object the conditions object
     * @return the object converted to a string
     */
    protected static String defaultToString(final BaseConditionsObject object) {
        final StringBuffer sb = new StringBuffer();
        sb.append(object.getClass().getSimpleName() + " { ");
        sb.append("id: " + object.getRowId() + ", ");
        for (final String field : object.getFieldValues().getFieldNames()) {
            sb.append(field + ": " + object.getFieldValue(Object.class, field) + ", ");
        }
        sb.setLength(sb.length() - 2);
        sb.append(" }");
        return sb.toString();
    }

    /**
     * The JDBC connection used for database operations.
     */
    private Connection connection;

    /**
     * The field values which is a map of key-value pairs corresponding to column values in the database.
     */
    private FieldValues fieldValues;

    /**
     * The row ID of the object in its table. This will be -1 for new objects that are not in the database.
     */
    private int rowId = UNSET_ROW_ID;

    /**
     * Flag to indicate that the object is locally changed and a database update has not been executed.
     */
    private boolean isDirty;

    /**
     * The information about the associated table such as the table and column names.
     */
    private TableMetaData tableMetaData;

    /**
     * No argument class constructor; usable by sub-classes only.
     */
    protected BaseConditionsObject() {
        this.fieldValues = new FieldValuesMap();
    }

    /**
     * Public class constructor.
     * <p>
     * This should be used when creating new objects without a list of field values. A new <code>FieldValues</code>
     * object will be automatically created from the table information.
     *
     * @param connection the database connection
     * @param tableMetaData the table meta data
     */
    public BaseConditionsObject(final Connection connection, final TableMetaData tableMetaData) {
        this.connection = connection;
        this.tableMetaData = tableMetaData;
        this.fieldValues = new FieldValuesMap(tableMetaData);
    }

    /**
     * Full qualified class constructor.
     * <p>
     * This should be used when creating new objects from a list of field values.
     *
     * @param connection the database connection
     * @param tableMetaData the table meta data
     * @param fields the field values
     */
    public BaseConditionsObject(final Connection connection, final TableMetaData tableMetaData, final FieldValues fields) {
        this.connection = connection;
        this.tableMetaData = tableMetaData;
        this.fieldValues = fields;

        // Since a list of field values are being provided, this object is flagged as needing to be updated.
        this.isDirty = true;
    }

    /**
     * Delete the object from the database using its row ID.
     * 
     * @throws DatabaseObjectException if object is not in the database
     * @throws SQLException if there is an error performing the delete operation
     */
    @Override
    public final void delete() throws DatabaseObjectException, SQLException {
        if (isNew()) {            
            throw new DatabaseObjectException("Object is not in the database.", this);
        }
        if (connection.getAutoCommit() == false) {
            connection.setAutoCommit(true);
        }
        PreparedStatement statement = null;
        try {
            statement = this.connection.prepareStatement("DELETE FROM " + this.tableMetaData.getTableName() + " WHERE id = ?");
            statement.setInt(1, + this.getRowId());
            statement.executeUpdate();
            this.rowId = UNSET_ROW_ID;
        } finally {
            if (statement != null) {
                statement.close();
            }
        }
    }

    /**
     * Get the collection ID of the object.
     * 
     * @return the collection ID of the object
     */
    @Override
    @Field(names = {"collection_id"})
    public final Integer getCollectionId() {
        if (this.fieldValues.isNonNull(COLLECTION_ID_FIELD)) {
            return getFieldValue(Integer.class, COLLECTION_ID_FIELD);
        } else {
            return UNSET_COLLECTION_ID;
        }
    }

    /**
     * Get the field values.
     * 
     * @return the field values
     */
    @Override
    public FieldValues getFieldValues() {
        return this.fieldValues;
    }

    /**
     * Get the row ID.
     * 
     * @return the row ID
     */
    @Override
    public final int getRowId() {
        return this.rowId;
    }

    /**
     * Get the table meta data for the object.
     * 
     * @return the table meta data or <code>null</code> if not set
     */
    @Override
    public final TableMetaData getTableMetaData() {
        return this.tableMetaData;
    }

    /**
     * Get a field value by name.
     * 
     * @param type the return type
     * @param name the name of the field
     */
    @Override
    public final <T> T getFieldValue(final Class<T> type, final String name) {
        return type.cast(this.fieldValues.getValue(type, name));
    }
   
    /**
     *
     */
    @Override
    // FIXME: Rewrite to use a PreparedStatement.
    public final void insert() throws DatabaseObjectException, SQLException {
        if (!this.isNew()) {
            throw new DatabaseObjectException("Cannot insert existing record with row ID: " + this.getRowId(), this);
        }
        if (!this.hasValidCollection()) {
            throw new DatabaseObjectException("Cannot insert object without a valid collection ID.", this);
        }
        final StringBuffer sb = new StringBuffer();
        sb.append("INSERT INTO " + this.tableMetaData.getTableName() + " (");
        for (final String fieldName : this.fieldValues.getFieldNames()) {
            sb.append(fieldName + ", ");
        }
        sb.setLength(sb.length() - 2);
        sb.append(") VALUES (");
        for (final Object value : this.fieldValues.getValues()) {
            if (value instanceof Date) {
                sb.append("STR_TO_DATE( '" + DATE_FORMAT.format((Date) value) + "', '%Y-%m-%d %H:%i:%S' ), ");
            } else {
                sb.append("'" + value + "', ");
            }
        }
        sb.setLength(sb.length() - 2);
        sb.append(")");
        final String sql = sb.toString();
        Statement statement = null;
        ResultSet resultSet = null;
        try {
            statement = this.connection.createStatement();
            statement.executeUpdate(sql, Statement.RETURN_GENERATED_KEYS);
            resultSet = statement.getGeneratedKeys();
            while (resultSet.next()) {
                final int key = resultSet.getInt(1);
                this.rowId = key;
                break;
            }
        } finally {
            if (resultSet != null) {
                resultSet.close();
            }
            if (statement != null) {
                statement.close();
            }
        }
        this.isDirty = false;
    }

    /**
     *
     */
    @Override
    public final boolean isDirty() {
        return this.isDirty;
    }

    /**
     *
     */
    @Override
    public final boolean isNew() {
        return getRowId() == UNSET_ROW_ID;
    }

    /**
     *
     */
    @Override
    public final boolean select(final int id) throws DatabaseObjectException, SQLException {
        this.rowId = id;
        if (id < 1) {
            throw new IllegalArgumentException("bad row ID value: " + id);
        }
        final StringBuffer sb = new StringBuffer();
        sb.append("SELECT");
        for (final String fieldName : this.tableMetaData.getFieldNames()) {
            sb.append(" " + fieldName + ",");
        }
        sb.setLength(sb.length() - 1);
        sb.append(" FROM " + this.tableMetaData.getTableName());
        sb.append(" WHERE id = " + this.getRowId());
        final String sql = sb.toString();
        Statement statement = null;
        ResultSet resultSet = null;
        boolean selected = false;
        try {
            statement = this.connection.createStatement();
            resultSet = statement.executeQuery(sql);
            selected = resultSet.next();
            if (selected) {
                for (int columnIndex = 1; columnIndex <= this.tableMetaData.getFieldNames().length; columnIndex++) {
                    this.setFieldValue(this.tableMetaData.getFieldNames()[columnIndex - 1],
                            resultSet.getObject(columnIndex));
                }
            }
        } finally {
            if (resultSet != null) {
                resultSet.close();
            }
            if (statement != null) {
                statement.close();
            }
        }
        return selected;
    }

    /**
     *
     */
    void setCollectionId(final int collectionId) throws ConditionsObjectException {
        //if (this.getCollectionId() != UNSET_COLLECTION_ID) {
        //    throw new ConditionsObjectException("The collection ID is already set on this object.");
        //}
        //if (collectionId <= UNSET_COLLECTION_ID) {
        //    throw new ConditionsObjectException("Invalid collection ID value: " + collectionId);
        //}
        this.setFieldValue(COLLECTION_ID_FIELD, collectionId);
    }

    /**
     *
     */
    @Override
    public final void setConnection(final Connection connection) {
        this.connection = connection;
    }

    /**
     *
     */
    @Override
    public void setFieldValues(final FieldValues fieldValues) {
        this.fieldValues = fieldValues;
        this.isDirty = true;
    }

    /**
     *
     */
    @Override
    public final void setTableMetaData(final TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }

    /**
     *
     */
    @Override
    public final void setFieldValue(final String name, final Object value) {
        this.fieldValues.setValue(name, value);
        this.isDirty = true;
    }

    /**
     *
     */
    @Override
    public String toString() {
        return defaultToString(this);
    }

    /**
     *
     */
    @Override
    public final boolean update() throws DatabaseObjectException, SQLException {
        int rowsUpdated = 0;
        if (isDirty()) {
            if (isNew()) {
                throw new DatabaseObjectException("Cannot update a new object.", this);
            }
            final StringBuffer sb = new StringBuffer();
            sb.append("UPDATE " + this.tableMetaData.getTableName() + " SET ");
            for (final String fieldName : this.fieldValues.getFieldNames()) {
                sb.append(fieldName + "=");
                final Object value = this.fieldValues.getValue(fieldName);
                if (value instanceof Date) {
                    // FIXME: Is there a more generic way to handle this?
                    sb.append("STR_TO_DATE( '" + DATE_FORMAT.format((Date) value) + "', '%Y-%m-%d %H:%i:%S' ), ");
                } else {
                    sb.append("'" + value + "', ");
                }
            }
            sb.setLength(sb.length() - 2);
            sb.append(" WHERE id=" + this.getRowId());
            final String sql = sb.toString();
            Statement statement = null;
            try {
                statement = this.connection.createStatement();
                rowsUpdated = statement.executeUpdate(sql);
            } finally {
                if (statement != null) {
                    statement.close();
                }
            }
        }
        if (rowsUpdated > 0) {
            this.isDirty = false;
        }
        return rowsUpdated != 0;
    }

    /**
     *
     */
    @Override
    public <T> T getFieldValue(final String name) {
        return (T) this.fieldValues.getValue(name);
    }

    /**
     *
     */
    @Override
    public boolean hasValidCollection() {
        return getCollectionId() != UNSET_COLLECTION_ID;
    }

    void setRowId(final int rowId) {
        this.rowId = rowId;
    }
}
