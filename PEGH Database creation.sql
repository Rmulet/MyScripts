CREATE DATABASE PEGH;
USE PEGH;

# Create the main tables of the database
CREATE table INTRAINDIVIDUAL (Tissue VARCHAR(30), Assembly VARCHAR(10), ID VARCHAR(20), Sex CHAR(1), Source VARCHAR(10));
CREATE table Interindividual (ExpId VARCHAR (10), Tissue VARCHAR(15), Individuals VARCHAR(5), Assembly VARCHAR(10), Source VARCHAR(15));
CREATE table Interspecies (ExpId VARCHAR (10), Tissue VARCHAR(15), Species VARCHAR(20), Individuals VARCHAR(5), Assemblies VARCHAR(10), Source VARCHAR(15));

# Some adjustments after table creation:
ALTER TABLE INTRAINDIVIDUAL RENAME TO Intraindividual;
ALTER TABLE Intraindividual 
CHANGE COLUMN ID
ID VARCHAR(20)
FIRST;
ALTER TABLE Intraindividual 
CHANGE COLUMN Tissue
Tissues VARCHAR(40);
ALTER TABLE Interindividual # Define the primary keyC
ADD PRIMARY KEY (ExpId); 
ALTER TABLE Intraindividual
ADD COLUMN ExpId VARCHAR(10)
FIRST; # AFTER S; 

# Create an example of data table. Iterind = Interindividual
CREATE table iterind1 (HS_Start VARCHAR (10), HS_End VARCHAR (10), MM_Start VARCHAR(10), MM_End VARCHAR(10), S VARCHAR(10), Theta VARCHAR(10), Pi VARCHAR (10), D VARCHAR(5), K VARCHAR(5), MKT VARCHAR(5), H3K27me3 VARCHAR(8),H3K4me3 VARCHAR(8), Cm VARCHAR(8));

# Adding text files to MySQL tables (mock data):
CREATE table trial1 (Chr CHAR(3), Source CHAR(10), Feature CHAR(10), Start VARCHAR(1), End VARCHAR(15), Score VARCHAR(10), Strand CHAR(1), Frame INT, Attribute VARCHAR(100));
LOAD DATA LOCAL INFILE '/home/roger/Documents/Output.txt' INTO TABLE trial1 COLUMNS TERMINATED BY '\t';

select * from iterind1;

# Table for genomics data:
CREATE table Genomics (Window VARCHAR(40), S VARCHAR(10), Pi VARCHAR(10), SFS VARCHAR(10), Divsites VARCHAR(5), P_ns VARCHAR(3), P_syn VARCHAR(3), D_ns VARCHAR(3), D_syn VARCHAR(3), NI VARCHAR(10), Unkown VARCHAR (3) );

# Add new rows:

INSERT INTO Genomics
VALUES

# View unique rows:

SELECT DISTINCT *
FROM Genomics
ORDER BY Window; 

# Add primary key to an existing table:

ALTER TABLE Genomics
CHANGE COLUMN Window
Window VARCHAR(20);

ALTER TABLE Genomics
ADD PRIMARY KEY (Window);SELECT * FROM PEGH.IHECMeta;

# Remove duplicates:

ALTER IGNORE TABLE jobs
ADD UNIQUE INDEX idx_name (site_id, title, company);

# Load a CSV table: (careful with the Primary Key!)

LOAD DATA LOCAL INFILE '/home/roger/Documents/3_EpigenomicsData/IHEC_Metadata.csv'
INTO TABLE IHECMeta
FIELDS TERMINATED BY ','
	ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 LINES;

http://www.mysqltutorial.org/import-csv-file-mysql-table/

# Delete records of a table:

TRUNCATE IHECMeta; # Or...
DELETE FROM PEGH.IHECMeta;

# Change one cell:

SET SQL_SAFE_UPDATES = 0;
UPDATE mytable
    SET column1 = value1,
        column2 = value2
    WHERE key_value = some_value;
SET SQL_SAFE_UPDATES = 1;

# "Escape" variable names: Enclose in `XX`



http://stackoverflow.com/questions/3024546/change-one-cells-data-in-mysql

