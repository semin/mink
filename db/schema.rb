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

ActiveRecord::Schema.define(:version => 0) do

  create_table :delayed_jobs, :force => true do |t|
    t.integer  :priority, :default => 0
    t.integer  :attempts, :default => 0
    t.text     :handler
    t.string   :last_error
    t.datetime :run_at
    t.datetime :locked_at
    t.datetime :failed_at
    t.string   :locked_by
    t.timestamps
  end


  # 'scop' table
  create_table :scop, :force => true do |t|
    t.belongs_to  :parent
    t.integer     :lft
    t.integer     :rgt
    t.string      :type
    t.integer     :sunid
    t.string      :stype
    t.string      :sccs
    t.string      :sid
    t.string      :description
    t.boolean     :rep95, :default => false
  end

  add_index :scop, :sunid
  add_index :scop, :parent_id
  add_index :scop, :lft
  add_index :scop, :rgt
  add_index :scop, [:id, :type]


  # 'mink_vectors' table
  create_table :mink_vectors, :force => true do |t|
    t.belongs_to  :scop
    t.string      :sid
    t.integer     :sunid
    t.string      :sccs
    t.float       :area_a
    t.float       :r_half_a
    t.float       :std_a
    t.float       :area_p
    t.float       :r_half_p
    t.float       :std_p
    t.float       :mean
    t.float       :std_mb
    t.float       :kurtosis
    t.float       :skewness
    t.float       :area_e
    t.float       :std_e
    t.float       :is
    t.string      :scop_class_description
    t.string      :scop_fold_description
    t.string      :scop_superfamily_description
    t.string      :scop_family_description
    t.string      :scop_protein_description
    t.string      :scop_species_description
    t.string      :scop_domain_description
  end

  add_index :mink_vectors, :sid,    :unique => true
  add_index :mink_vectors, :sunid,  :unique => true


  # 'norm_mink_vectors' table
  create_table :norm_mink_vectors, :force => true do |t|
    t.belongs_to  :scop
    t.belongs_to  :mink_vector
    t.string      :sid
    t.integer     :sunid
    t.float       :area_a
    t.float       :r_half_a
    t.float       :std_a
    t.float       :area_p
    t.float       :r_half_p
    t.float       :std_p
    t.float       :mean
    t.float       :std_mb
    t.float       :kurtosis
    t.float       :skewness
    t.float       :area_e
    t.float       :std_e
    t.float       :is
  end

  add_index :norm_mink_vectors, :sid,    :unique => true
  add_index :norm_mink_vectors, :sunid,  :unique => true


  # 'norm_mink_vectors' table
  create_table :mink_values, :force => true do |t|
    t.float       :min_area_a
    t.float       :min_r_half_a
    t.float       :min_std_a
    t.float       :min_area_p
    t.float       :min_r_half_p
    t.float       :min_std_p
    t.float       :min_mean
    t.float       :min_std_mb
    t.float       :min_kurtosis
    t.float       :min_skewness
    t.float       :min_area_e
    t.float       :min_std_e
    t.float       :min_is

    t.float       :max_area_a
    t.float       :max_r_half_a
    t.float       :max_std_a
    t.float       :max_area_p
    t.float       :max_r_half_p
    t.float       :max_std_p
    t.float       :max_mean
    t.float       :max_std_mb
    t.float       :max_kurtosis
    t.float       :max_skewness
    t.float       :max_area_e
    t.float       :max_std_e
    t.float       :max_is

    t.float       :submax_area_a
    t.float       :submax_r_half_a
    t.float       :submax_std_a
    t.float       :submax_area_p
    t.float       :submax_r_half_p
    t.float       :submax_std_p
    t.float       :submax_mean
    t.float       :submax_std_mb
    t.float       :submax_kurtosis
    t.float       :submax_skewness
    t.float       :submax_area_e
    t.float       :submax_std_e
    t.float       :submax_is
  end


  create_table :mink_vector_similarities, :force => true do |t|
    t.belongs_to  :mink_vector
    t.belongs_to  :similar_mink_vector
    t.float       :distance
  end

  add_index :mink_vector_similarities, :distance
  add_index :mink_vector_similarities, [:mink_vector_id, :similar_mink_vector_id], :name => "mink1_mink2"
  add_index :mink_vector_similarities, [:similar_mink_vector_id, :mink_vector_id], :name => "mink2_mink1"


  create_table :norm_mink_vector_similarities, :force => true do |t|
    t.belongs_to  :norm_mink_vector
    t.belongs_to  :similar_norm_mink_vector
    t.float       :distance
  end

  add_index :norm_mink_vector_similarities, :distance
  add_index :norm_mink_vector_similarities, [:norm_mink_vector_id, :similar_norm_mink_vector_id], :name => "norm_mink1_norm_mink2"
  add_index :norm_mink_vector_similarities, [:similar_norm_mink_vector_id, :norm_mink_vector_id], :name => "norm_mink2_norm_mink1"

  create_table :mink_searches, :force => true do |t|
    t.float       :cutoff
    t.string      :uuid
    t.string      :email
    t.timestamp   :started_at
    t.timestamp   :finished_at
    t.float       :elapsed_time
    t.string      :status
    t.float       :progress,     :default => 0
    t.float       :area_a
    t.float       :r_half_a
    t.float       :std_a
    t.float       :area_p
    t.float       :r_half_p
    t.float       :std_p
    t.float       :mean
    t.float       :std_mb
    t.float       :kurtosis
    t.float       :skewness
    t.float       :area_e
    t.float       :std_e
    t.float       :is
    t.float       :ori_area_a
    t.float       :ori_r_half_a
    t.float       :ori_std_a
    t.float       :ori_area_p
    t.float       :ori_r_half_p
    t.float       :ori_std_p
    t.float       :ori_mean
    t.float       :ori_std_mb
    t.float       :ori_kurtosis
    t.float       :ori_skewness
    t.float       :ori_area_e
    t.float       :ori_std_e
    t.float       :ori_is
    # columns for Paperclip plugin
    t.string      :pdb_file_name
    t.string      :pdb_content_type
    t.string      :pdb_file_size
    t.datetime    :pdb_updated_at
  end


  create_table :mink_search_hits, :force => true do |t|
    t.belongs_to  :mink_search
    t.belongs_to  :norm_mink_vector
    t.float       :distance
  end

  add_index :mink_search_hits, :distance


  create_table :news, :force => true do |t|
    t.date    :date
    t.string  :title
    t.text    :content
  end

  add_index :news, :date

end
