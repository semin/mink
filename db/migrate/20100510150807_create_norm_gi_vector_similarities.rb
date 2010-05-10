class CreateNormGiVectorSimilarities < ActiveRecord::Migration
  def self.up
    create_table :norm_gi_vector_similarities, :force => true do |t|
      t.belongs_to  :norm_gi_vector
      t.belongs_to  :similar_norm_gi_vector
      t.float       :distance
    end

    add_index :norm_gi_vector_similarities, :distance
    add_index :norm_gi_vector_similarities, [:norm_gi_vector_id, :similar_norm_gi_vector_id], :name => "norm_gi1_norm_gi2"
    add_index :norm_gi_vector_similarities, [:similar_norm_gi_vector_id, :norm_gi_vector_id], :name => "norm_gi2_norm_gi1"
  end

  def self.down
    drop_table :norm_gi_vector_similarities
  end
end
