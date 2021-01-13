def asedb_to_df(dbname,data={'id':'','formula':'','energy':'','natoms':'','volume':''}):
    import pandas as pd
    from ase.db import connect
    
    # connect to the database
    db = connect(dbname)
    # create an empty dataframe 
    df = pd.DataFrame()
    # add all rows into the dataframe
    for row in db.select():
        # extract only what is useful from default keys
        for col in data:
            data[col] = row.get(col)
        # extract all key-value pairs and add to the dictionary
        data.update(row.key_value_pairs)
        # append the dictionary into the dataframe
        df = df.append(pd.DataFrame([data]),ignore_index=True)
    return df
