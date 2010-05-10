# This file is auto-generated from the current state of the database. Instead of editing this file, 
# please use the migrations feature of Active Record to incrementally modify your database, and
# then regenerate this schema definition.
#
# Note that this schema.rb definition is the authoritative source for your database schema. If you need
# to create the application database on another system, you should be using db:schema:load, not running
# all the migrations from scratch. The latter is a flawed and unsustainable approach (the more migrations
# you'll amass, the slower it'll run and the greater likelihood for issues).
#
# It's strongly recommended to check this file into your version control system.

ActiveRecord::Schema.define(:version => 20100510150909) do

  create_table "delayed_jobs", :force => true do |t|
    t.integer  "priority",   :default => 0
    t.integer  "attempts",   :default => 0
    t.text     "handler"
    t.string   "last_error"
    t.datetime "run_at"
    t.datetime "locked_at"
    t.datetime "failed_at"
    t.string   "locked_by"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  create_table "gi_values", :force => true do |t|
    t.float "min_length"
    t.float "min_int12"
    t.float "min_inta12"
    t.float "min_int12_34"
    t.float "min_inta12_34"
    t.float "min_int12_a34"
    t.float "min_inta12_a34"
    t.float "min_int13_24"
    t.float "min_inta13_24"
    t.float "min_int13_a24"
    t.float "min_inta13_a24"
    t.float "min_int14_23"
    t.float "min_inta14_23"
    t.float "min_int14_a23"
    t.float "min_inta14_a23"
    t.float "min_int12_34_56"
    t.float "min_int12_35_46"
    t.float "min_int12_36_45"
    t.float "min_int13_24_56"
    t.float "min_int13_25_46"
    t.float "min_int13_26_45"
    t.float "min_int14_23_56"
    t.float "min_int14_25_36"
    t.float "min_int14_26_35"
    t.float "min_int15_23_46"
    t.float "min_int15_24_36"
    t.float "min_int15_26_34"
    t.float "min_int16_23_45"
    t.float "min_int16_24_35"
    t.float "min_int16_25_34"
    t.float "max_length"
    t.float "max_int12"
    t.float "max_inta12"
    t.float "max_int12_34"
    t.float "max_inta12_34"
    t.float "max_int12_a34"
    t.float "max_inta12_a34"
    t.float "max_int13_24"
    t.float "max_inta13_24"
    t.float "max_int13_a24"
    t.float "max_inta13_a24"
    t.float "max_int14_23"
    t.float "max_inta14_23"
    t.float "max_int14_a23"
    t.float "max_inta14_a23"
    t.float "max_int12_34_56"
    t.float "max_int12_35_46"
    t.float "max_int12_36_45"
    t.float "max_int13_24_56"
    t.float "max_int13_25_46"
    t.float "max_int13_26_45"
    t.float "max_int14_23_56"
    t.float "max_int14_25_36"
    t.float "max_int14_26_35"
    t.float "max_int15_23_46"
    t.float "max_int15_24_36"
    t.float "max_int15_26_34"
    t.float "max_int16_23_45"
    t.float "max_int16_24_35"
    t.float "max_int16_25_34"
    t.float "submax_length"
    t.float "submax_int12"
    t.float "submax_inta12"
    t.float "submax_int12_34"
    t.float "submax_inta12_34"
    t.float "submax_int12_a34"
    t.float "submax_inta12_a34"
    t.float "submax_int13_24"
    t.float "submax_inta13_24"
    t.float "submax_int13_a24"
    t.float "submax_inta13_a24"
    t.float "submax_int14_23"
    t.float "submax_inta14_23"
    t.float "submax_int14_a23"
    t.float "submax_inta14_a23"
    t.float "submax_int12_34_56"
    t.float "submax_int12_35_46"
    t.float "submax_int12_36_45"
    t.float "submax_int13_24_56"
    t.float "submax_int13_25_46"
    t.float "submax_int13_26_45"
    t.float "submax_int14_23_56"
    t.float "submax_int14_25_36"
    t.float "submax_int14_26_35"
    t.float "submax_int15_23_46"
    t.float "submax_int15_24_36"
    t.float "submax_int15_26_34"
    t.float "submax_int16_23_45"
    t.float "submax_int16_24_35"
    t.float "submax_int16_25_34"
  end

  create_table "gi_vector_similarities", :force => true do |t|
    t.integer "gi_vector_id"
    t.integer "similar_gi_vector_id"
    t.float   "distance"
  end

  add_index "gi_vector_similarities", ["distance"], :name => "index_gi_vector_similarities_on_distance"
  add_index "gi_vector_similarities", ["gi_vector_id", "similar_gi_vector_id"], :name => "gi1_gi2"
  add_index "gi_vector_similarities", ["similar_gi_vector_id", "gi_vector_id"], :name => "gi2_gi1"

  create_table "gi_vectors", :force => true do |t|
    t.integer "scop_id"
    t.string  "sid"
    t.integer "sunid"
    t.string  "sccs"
    t.string  "chain_code"
    t.integer "cas"
    t.integer "cas_missing"
    t.float   "length"
    t.float   "int12"
    t.float   "inta12"
    t.float   "int12_34"
    t.float   "inta12_34"
    t.float   "int12_a34"
    t.float   "inta12_a34"
    t.float   "int13_24"
    t.float   "inta13_24"
    t.float   "int13_a24"
    t.float   "inta13_a24"
    t.float   "int14_23"
    t.float   "inta14_23"
    t.float   "int14_a23"
    t.float   "inta14_a23"
    t.float   "int12_34_56"
    t.float   "int12_35_46"
    t.float   "int12_36_45"
    t.float   "int13_24_56"
    t.float   "int13_25_46"
    t.float   "int13_26_45"
    t.float   "int14_23_56"
    t.float   "int14_25_36"
    t.float   "int14_26_35"
    t.float   "int15_23_46"
    t.float   "int15_24_36"
    t.float   "int15_26_34"
    t.float   "int16_23_45"
    t.float   "int16_24_35"
    t.float   "int16_25_34"
    t.string  "scop_class_description"
    t.string  "scop_fold_description"
    t.string  "scop_superfamily_description"
    t.string  "scop_family_description"
    t.string  "scop_protein_description"
    t.string  "scop_species_description"
    t.string  "scop_domain_description"
  end

  add_index "gi_vectors", ["sccs"], :name => "class"
  add_index "gi_vectors", ["sid", "chain_code"], :name => "index_gi_vectors_on_sid_and_chain_code", :unique => true
  add_index "gi_vectors", ["sunid", "chain_code"], :name => "index_gi_vectors_on_sunid_and_chain_code", :unique => true

  create_table "gi_vectors_main4", :id => false, :force => true do |t|
    t.integer "id",                           :default => 0, :null => false
    t.integer "scop_id"
    t.string  "sid"
    t.integer "sunid"
    t.string  "sccs"
    t.string  "chain_code"
    t.integer "cas"
    t.integer "cas_missing"
    t.float   "length"
    t.float   "int12"
    t.float   "inta12"
    t.float   "int12_34"
    t.float   "inta12_34"
    t.float   "int12_a34"
    t.float   "inta12_a34"
    t.float   "int13_24"
    t.float   "inta13_24"
    t.float   "int13_a24"
    t.float   "inta13_a24"
    t.float   "int14_23"
    t.float   "inta14_23"
    t.float   "int14_a23"
    t.float   "inta14_a23"
    t.float   "int12_34_56"
    t.float   "int12_35_46"
    t.float   "int12_36_45"
    t.float   "int13_24_56"
    t.float   "int13_25_46"
    t.float   "int13_26_45"
    t.float   "int14_23_56"
    t.float   "int14_25_36"
    t.float   "int14_26_35"
    t.float   "int15_23_46"
    t.float   "int15_24_36"
    t.float   "int15_26_34"
    t.float   "int16_23_45"
    t.float   "int16_24_35"
    t.float   "int16_25_34"
    t.string  "scop_class_description"
    t.string  "scop_fold_description"
    t.string  "scop_superfamily_description"
    t.string  "scop_family_description"
    t.string  "scop_protein_description"
    t.string  "scop_species_description"
    t.string  "scop_domain_description"
  end

  create_table "git_values", :force => true do |t|
    t.float "min_measure1"
    t.float "min_measure2"
    t.float "min_measure3"
    t.float "min_measure4"
    t.float "min_measure5"
    t.float "min_measure6"
    t.float "min_measure7"
    t.float "min_measure8"
    t.float "min_measure9"
    t.float "min_measure10"
    t.float "min_measure11"
    t.float "min_measure12"
    t.float "min_measure13"
    t.float "min_measure14"
    t.float "min_measure15"
    t.float "min_measure16"
    t.float "min_measure17"
    t.float "min_measure18"
    t.float "min_measure19"
    t.float "min_measure20"
    t.float "min_measure21"
    t.float "min_measure22"
    t.float "min_measure23"
    t.float "min_measure24"
    t.float "min_measure25"
    t.float "min_measure26"
    t.float "min_measure27"
    t.float "min_measure28"
    t.float "min_measure29"
    t.float "min_measure30"
    t.float "max_measure1"
    t.float "max_measure2"
    t.float "max_measure3"
    t.float "max_measure4"
    t.float "max_measure5"
    t.float "max_measure6"
    t.float "max_measure7"
    t.float "max_measure8"
    t.float "max_measure9"
    t.float "max_measure10"
    t.float "max_measure11"
    t.float "max_measure12"
    t.float "max_measure13"
    t.float "max_measure14"
    t.float "max_measure15"
    t.float "max_measure16"
    t.float "max_measure17"
    t.float "max_measure18"
    t.float "max_measure19"
    t.float "max_measure20"
    t.float "max_measure21"
    t.float "max_measure22"
    t.float "max_measure23"
    t.float "max_measure24"
    t.float "max_measure25"
    t.float "max_measure26"
    t.float "max_measure27"
    t.float "max_measure28"
    t.float "max_measure29"
    t.float "max_measure30"
    t.float "submax_measure1"
    t.float "submax_measure2"
    t.float "submax_measure3"
    t.float "submax_measure4"
    t.float "submax_measure5"
    t.float "submax_measure6"
    t.float "submax_measure7"
    t.float "submax_measure8"
    t.float "submax_measure9"
    t.float "submax_measure10"
    t.float "submax_measure11"
    t.float "submax_measure12"
    t.float "submax_measure13"
    t.float "submax_measure14"
    t.float "submax_measure15"
    t.float "submax_measure16"
    t.float "submax_measure17"
    t.float "submax_measure18"
    t.float "submax_measure19"
    t.float "submax_measure20"
    t.float "submax_measure21"
    t.float "submax_measure22"
    t.float "submax_measure23"
    t.float "submax_measure24"
    t.float "submax_measure25"
    t.float "submax_measure26"
    t.float "submax_measure27"
    t.float "submax_measure28"
    t.float "submax_measure29"
    t.float "submax_measure30"
  end

  create_table "git_vector_similarities", :force => true do |t|
    t.integer "git_vector_id"
    t.integer "similar_git_vector_id"
    t.float   "distance"
  end

  add_index "git_vector_similarities", ["distance"], :name => "index_git_vector_similarities_on_distance"
  add_index "git_vector_similarities", ["git_vector_id", "similar_git_vector_id"], :name => "git1_git2"
  add_index "git_vector_similarities", ["similar_git_vector_id", "git_vector_id"], :name => "git2_git1"

  create_table "git_vectors", :force => true do |t|
    t.integer "scop_id"
    t.string  "sid"
    t.integer "sunid"
    t.string  "sccs"
    t.string  "chain_code"
    t.integer "cas_missing"
    t.integer "cas"
    t.float   "cube_root_cas_19_11"
    t.float   "measure1"
    t.float   "measure2"
    t.float   "measure3"
    t.float   "measure4"
    t.float   "measure5"
    t.float   "measure6"
    t.float   "measure7"
    t.float   "measure8"
    t.float   "measure9"
    t.float   "measure10"
    t.float   "measure11"
    t.float   "measure12"
    t.float   "measure13"
    t.float   "measure14"
    t.float   "measure15"
    t.float   "measure16"
    t.float   "measure17"
    t.float   "measure18"
    t.float   "measure19"
    t.float   "measure20"
    t.float   "measure21"
    t.float   "measure22"
    t.float   "measure23"
    t.float   "measure24"
    t.float   "measure25"
    t.float   "measure26"
    t.float   "measure27"
    t.float   "measure28"
    t.float   "measure29"
    t.float   "measure30"
    t.string  "scop_class_description"
    t.string  "scop_fold_description"
    t.string  "scop_superfamily_description"
    t.string  "scop_family_description"
    t.string  "scop_protein_description"
    t.string  "scop_species_description"
    t.string  "scop_domain_description"
  end

  add_index "git_vectors", ["sid", "chain_code"], :name => "index_git_vectors_on_sid_and_chain_code", :unique => true
  add_index "git_vectors", ["sunid", "chain_code"], :name => "index_git_vectors_on_sunid_and_chain_code", :unique => true

  create_table "mink_search_hits", :force => true do |t|
    t.integer "mink_search_id"
    t.integer "norm_mink_vector_id"
    t.float   "distance"
  end

  add_index "mink_search_hits", ["distance"], :name => "index_mink_search_hits_on_distance"

  create_table "mink_searches", :force => true do |t|
    t.float    "cutoff"
    t.string   "uuid"
    t.datetime "started_at"
    t.datetime "finished_at"
    t.float    "elapsed_time"
    t.string   "status"
    t.float    "progress",         :default => 0.0
    t.float    "area_a"
    t.float    "r_half_a"
    t.float    "std_a"
    t.float    "area_p"
    t.float    "r_half_p"
    t.float    "std_p"
    t.float    "mean"
    t.float    "std_mb"
    t.float    "kurtosis"
    t.float    "skewness"
    t.float    "area_e"
    t.float    "std_e"
    t.float    "is"
    t.float    "ori_area_a"
    t.float    "ori_r_half_a"
    t.float    "ori_std_a"
    t.float    "ori_area_p"
    t.float    "ori_r_half_p"
    t.float    "ori_std_p"
    t.float    "ori_mean"
    t.float    "ori_std_mb"
    t.float    "ori_kurtosis"
    t.float    "ori_skewness"
    t.float    "ori_area_e"
    t.float    "ori_std_e"
    t.float    "ori_is"
    t.string   "pdb_file_name"
    t.string   "pdb_content_type"
    t.string   "pdb_file_size"
    t.datetime "pdb_updated_at"
  end

  create_table "mink_values", :force => true do |t|
    t.float "min_area_a"
    t.float "min_r_half_a"
    t.float "min_std_a"
    t.float "min_area_p"
    t.float "min_r_half_p"
    t.float "min_std_p"
    t.float "min_mean"
    t.float "min_std_mb"
    t.float "min_kurtosis"
    t.float "min_skewness"
    t.float "min_area_e"
    t.float "min_std_e"
    t.float "min_is"
    t.float "max_area_a"
    t.float "max_r_half_a"
    t.float "max_std_a"
    t.float "max_area_p"
    t.float "max_r_half_p"
    t.float "max_std_p"
    t.float "max_mean"
    t.float "max_std_mb"
    t.float "max_kurtosis"
    t.float "max_skewness"
    t.float "max_area_e"
    t.float "max_std_e"
    t.float "max_is"
    t.float "submax_area_a"
    t.float "submax_r_half_a"
    t.float "submax_std_a"
    t.float "submax_area_p"
    t.float "submax_r_half_p"
    t.float "submax_std_p"
    t.float "submax_mean"
    t.float "submax_std_mb"
    t.float "submax_kurtosis"
    t.float "submax_skewness"
    t.float "submax_area_e"
    t.float "submax_std_e"
    t.float "submax_is"
  end

  create_table "mink_vector_similarities", :force => true do |t|
    t.integer "mink_vector_id"
    t.integer "similar_mink_vector_id"
    t.float   "distance"
  end

  add_index "mink_vector_similarities", ["distance"], :name => "index_mink_vector_similarities_on_distance"
  add_index "mink_vector_similarities", ["mink_vector_id", "similar_mink_vector_id"], :name => "mink1_mink2"
  add_index "mink_vector_similarities", ["similar_mink_vector_id", "mink_vector_id"], :name => "mink2_mink1"

  create_table "mink_vectors", :force => true do |t|
    t.integer "scop_id"
    t.string  "sid"
    t.integer "sunid"
    t.string  "sccs"
    t.float   "area_a"
    t.float   "r_half_a"
    t.float   "std_a"
    t.float   "area_p"
    t.float   "r_half_p"
    t.float   "std_p"
    t.float   "mean"
    t.float   "std_mb"
    t.float   "kurtosis"
    t.float   "skewness"
    t.float   "area_e"
    t.float   "std_e"
    t.float   "is"
    t.string  "scop_class_description"
    t.string  "scop_fold_description"
    t.string  "scop_superfamily_description"
    t.string  "scop_family_description"
    t.string  "scop_protein_description"
    t.string  "scop_species_description"
    t.string  "scop_domain_description"
  end

  add_index "mink_vectors", ["sccs"], :name => "class"
  add_index "mink_vectors", ["sid"], :name => "index_mink_vectors_on_sid", :unique => true
  add_index "mink_vectors", ["sunid"], :name => "index_mink_vectors_on_sunid", :unique => true

  create_table "mink_vectors_main4", :id => false, :force => true do |t|
    t.integer "id",                           :default => 0, :null => false
    t.integer "scop_id"
    t.string  "sid"
    t.integer "sunid"
    t.string  "sccs"
    t.float   "area_a"
    t.float   "r_half_a"
    t.float   "std_a"
    t.float   "area_p"
    t.float   "r_half_p"
    t.float   "std_p"
    t.float   "mean"
    t.float   "std_mb"
    t.float   "kurtosis"
    t.float   "skewness"
    t.float   "area_e"
    t.float   "std_e"
    t.float   "is"
    t.string  "scop_class_description"
    t.string  "scop_fold_description"
    t.string  "scop_superfamily_description"
    t.string  "scop_family_description"
    t.string  "scop_protein_description"
    t.string  "scop_species_description"
    t.string  "scop_domain_description"
  end

  create_table "news", :force => true do |t|
    t.date   "date"
    t.string "title"
    t.text   "content"
  end

  add_index "news", ["date"], :name => "index_news_on_date"

  create_table "norm_gi_vector_similarities", :force => true do |t|
    t.integer "norm_gi_vector_id"
    t.integer "similar_norm_gi_vector_id"
    t.float   "distance"
  end

  add_index "norm_gi_vector_similarities", ["distance"], :name => "index_norm_gi_vector_similarities_on_distance"
  add_index "norm_gi_vector_similarities", ["norm_gi_vector_id", "similar_norm_gi_vector_id"], :name => "norm_gi1_norm_gi2"
  add_index "norm_gi_vector_similarities", ["similar_norm_gi_vector_id", "norm_gi_vector_id"], :name => "norm_gi2_norm_gi1"

  create_table "norm_gi_vectors", :force => true do |t|
    t.integer "gi_vector_id"
    t.integer "scop_id"
    t.string  "sid"
    t.integer "sunid"
    t.string  "sccs"
    t.string  "chain_code"
    t.integer "cas"
    t.integer "cas_missing"
    t.float   "length"
    t.float   "int12"
    t.float   "inta12"
    t.float   "int12_34"
    t.float   "inta12_34"
    t.float   "int12_a34"
    t.float   "inta12_a34"
    t.float   "int13_24"
    t.float   "inta13_24"
    t.float   "int13_a24"
    t.float   "inta13_a24"
    t.float   "int14_23"
    t.float   "inta14_23"
    t.float   "int14_a23"
    t.float   "inta14_a23"
    t.float   "int12_34_56"
    t.float   "int12_35_46"
    t.float   "int12_36_45"
    t.float   "int13_24_56"
    t.float   "int13_25_46"
    t.float   "int13_26_45"
    t.float   "int14_23_56"
    t.float   "int14_25_36"
    t.float   "int14_26_35"
    t.float   "int15_23_46"
    t.float   "int15_24_36"
    t.float   "int15_26_34"
    t.float   "int16_23_45"
    t.float   "int16_24_35"
    t.float   "int16_25_34"
    t.string  "scop_class_description"
    t.string  "scop_fold_description"
    t.string  "scop_superfamily_description"
    t.string  "scop_family_description"
    t.string  "scop_protein_description"
    t.string  "scop_species_description"
    t.string  "scop_domain_description"
  end

  add_index "norm_gi_vectors", ["sid", "chain_code"], :name => "index_norm_gi_vectors_on_sid_and_chain_code", :unique => true
  add_index "norm_gi_vectors", ["sunid", "chain_code"], :name => "index_norm_gi_vectors_on_sunid_and_chain_code", :unique => true

  create_table "norm_git_vector_similarities", :force => true do |t|
    t.integer "norm_git_vector_id"
    t.integer "similar_norm_git_vector_id"
    t.float   "distance"
  end

  add_index "norm_git_vector_similarities", ["distance"], :name => "index_norm_git_vector_similarities_on_distance"
  add_index "norm_git_vector_similarities", ["norm_git_vector_id", "similar_norm_git_vector_id"], :name => "norm_git1_norm_git2"
  add_index "norm_git_vector_similarities", ["similar_norm_git_vector_id", "norm_git_vector_id"], :name => "norm_git2_norm_git1"

  create_table "norm_git_vectors", :force => true do |t|
    t.integer "git_vector_id"
    t.integer "scop_id"
    t.string  "sid"
    t.integer "sunid"
    t.string  "sccs"
    t.string  "chain_code"
    t.integer "cas_missing"
    t.integer "cas"
    t.float   "cube_root_cas_19_11"
    t.float   "measure1"
    t.float   "measure2"
    t.float   "measure3"
    t.float   "measure4"
    t.float   "measure5"
    t.float   "measure6"
    t.float   "measure7"
    t.float   "measure8"
    t.float   "measure9"
    t.float   "measure10"
    t.float   "measure11"
    t.float   "measure12"
    t.float   "measure13"
    t.float   "measure14"
    t.float   "measure15"
    t.float   "measure16"
    t.float   "measure17"
    t.float   "measure18"
    t.float   "measure19"
    t.float   "measure20"
    t.float   "measure21"
    t.float   "measure22"
    t.float   "measure23"
    t.float   "measure24"
    t.float   "measure25"
    t.float   "measure26"
    t.float   "measure27"
    t.float   "measure28"
    t.float   "measure29"
    t.float   "measure30"
    t.string  "scop_class_description"
    t.string  "scop_fold_description"
    t.string  "scop_superfamily_description"
    t.string  "scop_family_description"
    t.string  "scop_protein_description"
    t.string  "scop_species_description"
    t.string  "scop_domain_description"
  end

  add_index "norm_git_vectors", ["sid", "chain_code"], :name => "index_norm_git_vectors_on_sid_and_chain_code", :unique => true
  add_index "norm_git_vectors", ["sunid", "chain_code"], :name => "index_norm_git_vectors_on_sunid_and_chain_code", :unique => true

  create_table "norm_mink_vector_similarities", :force => true do |t|
    t.integer "norm_mink_vector_id"
    t.integer "similar_norm_mink_vector_id"
    t.float   "distance"
  end

  add_index "norm_mink_vector_similarities", ["distance"], :name => "index_norm_mink_vector_similarities_on_distance"
  add_index "norm_mink_vector_similarities", ["norm_mink_vector_id", "similar_norm_mink_vector_id"], :name => "norm_mink1_norm_mink2"
  add_index "norm_mink_vector_similarities", ["similar_norm_mink_vector_id", "norm_mink_vector_id"], :name => "norm_mink2_norm_mink1"

  create_table "norm_mink_vectors", :force => true do |t|
    t.integer "scop_id"
    t.integer "mink_vector_id"
    t.string  "sid"
    t.integer "sunid"
    t.string  "sccs"
    t.float   "area_a"
    t.float   "r_half_a"
    t.float   "std_a"
    t.float   "area_p"
    t.float   "r_half_p"
    t.float   "std_p"
    t.float   "mean"
    t.float   "std_mb"
    t.float   "kurtosis"
    t.float   "skewness"
    t.float   "area_e"
    t.float   "std_e"
    t.float   "is"
    t.string  "scop_class_description"
    t.string  "scop_fold_description"
    t.string  "scop_superfamily_description"
    t.string  "scop_family_description"
    t.string  "scop_protein_description"
    t.string  "scop_species_description"
    t.string  "scop_domain_description"
  end

  add_index "norm_mink_vectors", ["sid"], :name => "index_norm_mink_vectors_on_sid", :unique => true
  add_index "norm_mink_vectors", ["sunid"], :name => "index_norm_mink_vectors_on_sunid", :unique => true

  create_table "scop", :primary_key => "px", :force => true do |t|
    t.string  "sid"
    t.string  "pdbid",  :null => false
    t.string  "sccs",   :null => false
    t.integer "fam",    :null => false
    t.integer "supfam", :null => false
    t.integer "fold",   :null => false
    t.integer "class",  :null => false
    t.integer "res10"
  end

  add_index "scop", ["fam"], :name => "family"
  add_index "scop", ["fold"], :name => "class"
  add_index "scop", ["fold"], :name => "fold"
  add_index "scop", ["pdbid"], :name => "pdbid"
  add_index "scop", ["res10"], :name => "resolution"
  add_index "scop", ["supfam"], :name => "superfam"

  create_table "scops", :force => true do |t|
    t.integer "parent_id"
    t.integer "lft"
    t.integer "rgt"
    t.string  "type"
    t.integer "sunid"
    t.string  "stype"
    t.string  "sccs"
    t.string  "sid"
    t.string  "description"
    t.boolean "rep95",       :default => false
  end

  add_index "scops", ["id", "type"], :name => "index_scops_on_id_and_type"
  add_index "scops", ["lft"], :name => "index_scops_on_lft"
  add_index "scops", ["parent_id"], :name => "index_scops_on_parent_id"
  add_index "scops", ["rgt"], :name => "index_scops_on_rgt"
  add_index "scops", ["sunid"], :name => "index_scops_on_sunid"

end
