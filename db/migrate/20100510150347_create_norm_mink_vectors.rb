class CreateNormMinkVectors < ActiveRecord::Migration
  def self.up
  create_table :norm_mink_vectors, :force => true do |t|
    t.belongs_to  :scop
    t.belongs_to  :mink_vector
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

  add_index :norm_mink_vectors, :sid,    :unique => true
  add_index :norm_mink_vectors, :sunid,  :unique => true
  end

  def self.down
    drop_table :norm_mink_vectors
  end
end
