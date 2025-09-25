CREATE TABLE compound_candidate (
	id INTEGER NOT NULL, 
	confidence_rank INTEGER, 
	structure_per_id_rank INTEGER, 
	formula_rank INTEGER, 
	num_adducts INTEGER, 
	num_predicted_fingerprints INTEGER, 
	confidence_score FLOAT, 
	finger_id_score FLOAT, 
	zodiac_score FLOAT, 
	sirius_score FLOAT, 
	formula_sirius VARCHAR, 
	adduct_sirius VARCHAR, 
	inchi VARCHAR, 
	name_sirius VARCHAR, 
	smiles VARCHAR, 
	xlogp FLOAT, 
	rt_seconds FLOAT, 
	sirius_compound_folder VARCHAR, 
	feature_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(feature_id) REFERENCES feature_sirius (feature_id)
);


CREATE TABLE feature_sirius (
	feature_id INTEGER NOT NULL, 
	use_zodiac_scoring_for_best BOOLEAN NOT NULL, 
	highest_scoring_formula VARCHAR, 
	combined_feature_id INTEGER, 
	PRIMARY KEY (feature_id), 
	FOREIGN KEY(combined_feature_id) REFERENCES features (id)
);


CREATE TABLE features (
	id INTEGER NOT NULL, 
	feature_id INTEGER, 
	PRIMARY KEY (id)
);


CREATE TABLE compound_group (
	id INTEGER NOT NULL, 
	sirius_compound_folder VARCHAR, 
	formula_sirius VARCHAR, 
	adduct_sirius VARCHAR, 
	npc_pathway_name VARCHAR, 
	npc_pathway_probability FLOAT, 
	npc_superclass_name VARCHAR, 
	npc_superclass_probability FLOAT, 
	npc_class_name VARCHAR, 
	npc_class_probability FLOAT, 
	cf_most_specific_name VARCHAR, 
	cf_most_specific_probability FLOAT, 
	cf_level5_name VARCHAR, 
	cf_level5_probability FLOAT, 
	cf_subclass_name VARCHAR, 
	cf_subclass_probability FLOAT, 
	cf_class_name VARCHAR, 
	cf_class_probability FLOAT, 
	cf_superclass_name VARCHAR, 
	cf_superclass_probability FLOAT, 
	cf_path VARCHAR, 
	feature_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(feature_id) REFERENCES feature_sirius (feature_id)
);


CREATE TABLE formula_candidate (
	id INTEGER NOT NULL, 
	formula_sirius VARCHAR, 
	formula_rank INTEGER, 
	adduct_sirius VARCHAR, 
	zodiac_score FLOAT, 
	sirius_score FLOAT, 
	tree_score FLOAT, 
	isotope_score FLOAT, 
	num_explained_peaks INTEGER, 
	explained_intensity FLOAT, 
	median_mass_error_fragments_ppm FLOAT, 
	mass_error_precursor_ppm FLOAT, 
	lipid_class VARCHAR, 
	rt_seconds INTEGER, 
	sirius_compound_folder VARCHAR, 
	feature_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(feature_id) REFERENCES feature_sirius (feature_id)
);


CREATE TABLE gnps_features (
	id INTEGER NOT NULL, 
	feature_id INTEGER, 
	cluster_label INTEGER, 
	"M_gnps" FLOAT, 
	rt_seconds FLOAT, 
	other VARCHAR, 
	combined_feature_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(combined_feature_id) REFERENCES features (id)
);


CREATE TABLE intensities (
	id INTEGER NOT NULL, 
	sample_name VARCHAR, 
	value INTEGER, 
	feature_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(feature_id) REFERENCES metaboscape_features (id)
);


CREATE TABLE metaboscape_features (
	id INTEGER NOT NULL, 
	feature_id INTEGER, 
	rt_seconds FLOAT, 
	"M_metaboscape" FLOAT, 
	"CCS" FLOAT, 
	sigma_score FLOAT, 
	name_metaboscape VARCHAR, 
	formula_metaboscape VARCHAR, 
	adduct_metaboscape VARCHAR, 
	"KEGG" FLOAT, 
	"CAS" FLOAT, 
	combined_feature_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(combined_feature_id) REFERENCES features (id)
);


CREATE TABLE mgf_features (
	id INTEGER NOT NULL, 
	feature_id INTEGER, 
	polarity VARCHAR(3), 
	has_multiple_adducts BOOLEAN, 
	mz FLOAT, 
	charge INTEGER, 
	rt_seconds FLOAT, 
	ion VARCHAR, 
	rt_minutes FLOAT, 
	ms1_id INTEGER, 
	ms2_id INTEGER, 
	combined_feature_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(ms1_id) REFERENCES peak_list (id), 
	FOREIGN KEY(ms2_id) REFERENCES peak_list (id), 
	FOREIGN KEY(combined_feature_id) REFERENCES features (id)
);


CREATE TABLE peak_list (
	id INTEGER NOT NULL, 
	name VARCHAR, 
	PRIMARY KEY (id)
);


CREATE TABLE ms_spec (
	id INTEGER NOT NULL, 
	mz FLOAT, 
	ms_level INTEGER, 
	charge INTEGER, 
	rt_seconds FLOAT, 
	ion VARCHAR, 
	rt_minutes FLOAT, 
	peaks_id INTEGER, 
	feature_mgf_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(peaks_id) REFERENCES peak_list (id), 
	FOREIGN KEY(feature_mgf_id) REFERENCES mgf_features (id)
);


CREATE TABLE peak (
	id INTEGER NOT NULL, 
	mz FLOAT NOT NULL, 
	rt FLOAT, 
	intensity FLOAT, 
	peak_list_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(peak_list_id) REFERENCES peak_list (id)
);
