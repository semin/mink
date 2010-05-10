class CreateNormGitVectors < ActiveRecord::Migration
  def self.up
    create_table :norm_git_vectors, :force => true do |t|
      t.belongs_to  :git_vector
      t.belongs_to  :scop
      t.string      :sid
      t.integer     :sunid
      t.string      :sccs
      t.string      :chain_code
      t.integer     :cas_missing
      t.integer     :cas
      t.float       :cube_root_cas_19_11
      t.float       :measure1
      t.float       :measure2
      t.float       :measure3
      t.float       :measure4
      t.float       :measure5
      t.float       :measure6
      t.float       :measure7
      t.float       :measure8
      t.float       :measure9
      t.float       :measure10
      t.float       :measure11
      t.float       :measure12
      t.float       :measure13
      t.float       :measure14
      t.float       :measure15
      t.float       :measure16
      t.float       :measure17
      t.float       :measure18
      t.float       :measure19
      t.float       :measure20
      t.float       :measure21
      t.float       :measure22
      t.float       :measure23
      t.float       :measure24
      t.float       :measure25
      t.float       :measure26
      t.float       :measure27
      t.float       :measure28
      t.float       :measure29
      t.float       :measure30
      t.string      :scop_class_description
      t.string      :scop_fold_description
      t.string      :scop_superfamily_description
      t.string      :scop_family_description
      t.string      :scop_protein_description
      t.string      :scop_species_description
      t.string      :scop_domain_description
    end

    add_index :norm_git_vectors, [:sid, :chain_code],   :unique => true
    add_index :norm_git_vectors, [:sunid, :chain_code], :unique => true
  end

  def self.down
    drop_table :norm_git_vectors
  end
end
