class CreateNormMinkVectorSimilarities < ActiveRecord::Migration
  def self.up
    create_table :norm_mink_vector_similarities, :force => true do |t|
      t.belongs_to  :norm_mink_vector
      t.belongs_to  :similar_norm_mink_vector
      t.float       :distance
    end

    add_index :norm_mink_vector_similarities, :distance
    add_index :norm_mink_vector_similarities, [:norm_mink_vector_id, :similar_norm_mink_vector_id], :name => "norm_mink1_norm_mink2"
    add_index :norm_mink_vector_similarities, [:similar_norm_mink_vector_id, :norm_mink_vector_id], :name => "norm_mink2_norm_mink1"
  end

  def self.down
    drop_table :norm_mink_vector_similarities
  end
end
