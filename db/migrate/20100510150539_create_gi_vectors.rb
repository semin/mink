class CreateGiVectors < ActiveRecord::Migration
  def self.up
    create_table :gi_vectors, :force => true do |t|
      t.belongs_to  :scop
      t.string      :sid
      t.integer     :sunid
      t.string      :sccs
      t.string      :chain_code
      t.integer     :cas
      t.integer     :cas_missing
      t.float       :length
      t.float       :int12
      t.float       :inta12
      t.float       :int12_34
      t.float       :inta12_34
      t.float       :int12_a34
      t.float       :inta12_a34
      t.float       :int13_24
      t.float       :inta13_24
      t.float       :int13_a24
      t.float       :inta13_a24
      t.float       :int14_23
      t.float       :inta14_23
      t.float       :int14_a23
      t.float       :inta14_a23
      t.float       :int12_34_56
      t.float       :int12_35_46
      t.float       :int12_36_45
      t.float       :int13_24_56
      t.float       :int13_25_46
      t.float       :int13_26_45
      t.float       :int14_23_56
      t.float       :int14_25_36
      t.float       :int14_26_35
      t.float       :int15_23_46
      t.float       :int15_24_36
      t.float       :int15_26_34
      t.float       :int16_23_45
      t.float       :int16_24_35
      t.float       :int16_25_34
      t.string      :scop_class_description
      t.string      :scop_fold_description
      t.string      :scop_superfamily_description
      t.string      :scop_family_description
      t.string      :scop_protein_description
      t.string      :scop_species_description
      t.string      :scop_domain_description
    end

    add_index :gi_vectors, [:sid, :chain_code],   :unique => true
    add_index :gi_vectors, [:sunid, :chain_code], :unique => true
  end

  def self.down
    drop_table :gi_vectors
  end
end
