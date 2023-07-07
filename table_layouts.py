"""
All tables used in the initial database.
"""

compound_main = """ CREATE TABLE IF NOT EXISTS compound_main( 
            compound_id INTEGER PRIMARY KEY, 
            smiles TEXT NOT NULL,
            png TEXT NOT NULL, 
            volume REAL NOT NULL,
            concentration REAL NOT NULL,
            ac_id INTEGER NOT NULL,
            active TEXT NOT NULL,
            location INTEGER NOT NULL,
            FOREIGN KEY (ac_id) REFERENCES origin(ac_id),
            FOREIGN KEY (location) REFERENCES location(loc_id)
            ); """

compound_mp_table = """ CREATE TABLE IF NOT EXISTS compound_mp( 
            row_counter INTEGER PRIMARY KEY,
            mp_barcode TEXT NOT NULL UNIQUE,
            compound_id INTEGER NOT NULL,
            mp_well TEXT NOT NULL,
            volume REAL NOT NULL,
            date REAL NOT NULL,
            active INTEGER NOT NULL,
            freeze_thaw INTEGER NOT NULL,
            plate_type TEXT NOT NULL,
            location INTEGER NOT NULL,
            FOREIGN KEY (mp_barcode) REFERENCES mp_plates(mp_barcode),
            FOREIGN KEY (compound_id) REFERENCES compound_main(compound_id),
            FOREIGN KEY (location) REFERENCES location(loc_id),
            FOREIGN KEY (location) REFERENCES location(loc_id),
            FOREIGN KEY (plate_type) REFERENCES plate_type(plate_type)
            ); """

compound_dp_table = """ CREATE TABLE IF NOT EXISTS compound_dp( 
            row_counter INTEGER PRIMARY KEY,
            dp_barcode TEXT NOT NULL UNIQUE,
            compound_id INTEGER NOT NULL,
            dp_well TEXT NOT NULL, 
            volume REAL NOT NULL, 
            date REAL NOT NULL,         
            active TEXT NOT NULL,
            freeze_thaw INTEGER NOT NULL,
            plate_type TEXT NOT NULL,
            location INTEGER NOT NULL,
            mp_barcode TEXT NOT NULL,
            mp_well TEXT NOT NULL,
            FOREIGN KEY (compound_id) REFERENCES compound_main(compound_id),
            FOREIGN KEY (mp_barcode) REFERENCES mp_plates(mp_barcode),
            FOREIGN KEY (dp_barcode) REFERENCES dp_plates(dp_barcode),
            FOREIGN KEY (location) REFERENCES location(loc_id),
            FOREIGN KEY (plate_type) REFERENCES plate_type(plate_type)
            ); """

plate_type = """ CREATE TABLE IF NOT EXISTS plate_type( 
            row_id	INTEGER PRIMARY KEY,
            plate_type	TEXT NOT NULL UNIQUE,
            vendor	TEXT NOT NULL,
            product_number	TEXT NOT NULL,
            sterile	INTEGER NOT NULL,
            info	BLOB,
            size INTEGER
            well_offset_x	INTEGER,
            well_offset_y	INTEGER,
            well_spacing_x	INTEGER,
            well_spacing_y	INTEGER,
            plate_height	INTEGER,
            plate_height_lid	INTEGER,
            flang_height	INTEGER,
            well_depth	INTEGER,
            well_width	INTEGER,
            max_volume	INTEGER,
            working_volume	INTEGER,
            dead_volume	INTEGER,
            FOREIGN KEY("vendor") REFERENCES "vendors"("name")
            """

compound_data_table = """ CREATE TABLE IF NOT EXISTS compound_data( 
            compound_id INTEGER NOT NULL,
            exp_id INTEGER NOT NULL UNIQUE,
            FOREIGN KEY (compound_id) REFERENCES compound_main(compound_id),
            FOREIGN KEY (exp_id) REFERENCES experiment(exp_id)
            ); """

mother_plate_table = """ CREATE TABLE IF NOT EXISTS mp_plates(
            mp_barcode TEXT PRIMARY KEY,
            date REAL NOT NULL
            ); """

daughter_plate_table = """ CREATE TABLE IF NOT EXISTS dp_plates(
            dp_barcode TEXT PRIMARY KEY,
            date REAL NOT NULL
            ); """

location_table = """ CREATE TABLE IF NOT EXISTS locations(              
            loc_id INTEGER PRIMARY KEY AUTOINCREMENT,
            room TEXT NOT NULL,
            building TEXT NOT NULL UNIQUE,
            spot TEXT NOT NULL
            ); """

bio_experiment_table = """ CREATE TABLE IF NOT EXISTS bio_experiment(
            exp_id INTEGER PRIMARY KEY AUTOINCREMENT,
            assay_name TEXT NOT NULL,
            raw_data TEXT NOT NULL,
            responsible TEXT,
            date REAL NOT NULL,
            plate_name	TEXT NOT NULL UNIQUE,
            approval	TEXT NOT NULL,
            z_prime	REAL,
            process_data    TEXT,
            plate_layout	TEXT NOT NULL,
            note    TEXT,
            FOREIGN KEY("assay_name") REFERENCES "assay"("assay_name")
            ); """

biological_compound_data = """ CREATE TABLE IF NOT EXISTS biological_compound_data(
            bio_data_id INTEGER PRIMARY KEY UNIQUE,
            compound_id INTEGER NOT NULL,
            exp_id INTEGER NOT NULL,
            raw_data REAL NOT NULL,
            score REAL NOT NULL,
            hit TEXT NOT NULL,
            concentration   REAL NOT NULL,
            approved TEXT NOT NULL,
            note TEXT,
            FOREIGN KEY (compound_id) REFERENCES compound_main(compound_id),
            FOREIGN KEY (exp_id) REFERENCES exp_id(bio_experiment)
            ); """

lc_experiment_table = """ CREATE TABLE IF NOT EXISTS lc_experiment(
            row_id INTEGER PRIMARY KEY,
            batch TEXT NOT NULL,
            date REAL NOT NULL
            ); """

purity_data = """ CREATE TABLE IF NOT EXISTS purity(
            purity_id INTEGER PRIMARY KEY,
            compound_id INTEGER NOT NULL,
            batch TEXT NOT NULL,
            result_max REAL NOT NULL,
            result_max_ion REAL NOT NULL,
            result_total TEXT NOT NULL,
            date REAL NOT NULL,
            FOREIGN KEY (compound_id) REFERENCES compound_main(compound_id),
            FOREIGN KEY (batch) REFERENCES lc_experiment(batch)
            ); """

lcms_experiment_raw = """ CREATE TABLE IF NOT EXISTS lc_raw(
            row_id INTEGER PRIMARY KEY,
            compound_id INTEGER NOT NULL,
            sample TEXT NOT NULL,
            batch TEXT NOT NULL,
            method TEXT NOT NULL,
            file_name TEXT NOT NULL,
            date REAL NOT NULL
            FOREIGN KEY (compound_id) REFERENCES compound_main(compound_id),
            FOREIGN KEY (batch) REFERENCES lc_experiment(batch)
            ); """

compound_source = """ CREATE TABLE IF NOT EXISTS origin(
            ac_id INTEGER PRIMARY KEY AUTOINCREMENT,
            origin TEXT NOT NULL UNIQUE,
            ac TEXT NOT NULL,
            contact_person TEXT,
            e_mail TEXT,
            info BLOB
            ); """

customers = """ CREATE TABLE IF NOT EXISTS customers(
            row_id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL UNIQUE,
            e_mail TEXT NOT NULL,
            info BLOB
            ); """

vendors = """ CREATE TABLE IF NOT EXISTS vendors(
            row_id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL UNIQUE,
            e_mail REAL NOT NULL,
            info BLOB
            ); """

responsible = """ CREATE TABLE IF NOT EXISTS responsible(
            row_id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL UNIQUE,
            e_mail REAL NOT NULL,
            info BLOB
            ); """

assay = """ CREATE TABLE IF NOT EXISTS assay(
            row_id INTEGER PRIMARY KEY AUTOINCREMENT,         
            assay TEXT NOT NULL UNIQUE,
            plate_layout TEXT,
            sop TEXT,
            plates_run INTEGER,
            compounds_run INTEGER,
            z_prime_threshold	REAL,
            hit_threshold	REAL,
            FOREIGN KEY(plate_layout) REFERENCES plate_layout(plate_name)
            ); """

assay_customers = """ CREATE TABLE IF NOT EXISTS assay_customers(
            row_id INTEGER PRIMARY KEY AUTOINCREMENT,
            customer TEXT NOT NULL,
            exp_id INTEGER NOT NULL, 
            FOREIGN KEY (customer) REFERENCES customers(name),
            FOREIGN KEY (exp_id) REFERENCES bio_experiment(exp_id)
            ); """

plate_layout = """ CREATE TABLE IF NOT EXISTS plate_layout(  
            plate_name TEXT NOT NULL UNIQUE,
            well_layout TEXT,
            plate_type TEXT
            ); """