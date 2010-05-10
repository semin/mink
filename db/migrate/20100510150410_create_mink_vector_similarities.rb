class CreateMinkVectorSimilarities < ActiveRecord::Migration
  def self.up
    create_table :mink_vector_similarities, :force => true do |t|
      t.belongs_to  :mink_vector
      t.belongs_to  :similar_mink_vector
      t.float       :distance
    end

    add_index :mink_vector_similarities, :distance
    add_index :mink_vector_similarities, [:mink_vector_id, :similar_mink_vector_id], :name => "mink1_mink2"
    add_index :mink_vector_similarities, [:similar_mink_vector_id, :mink_vector_id], :name => "mink2_mink1"
  end

  def self.down
    drop_table mink_vector_similarities
  end
end
