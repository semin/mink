class CreateGiVectorSimilarities < ActiveRecord::Migration
  def self.up
    create_table :gi_vector_similarities, :force => true do |t|
      t.belongs_to  :gi_vector
      t.belongs_to  :similar_gi_vector
      t.float       :distance
    end

    add_index :gi_vector_similarities, :distance
    add_index :gi_vector_similarities, [:gi_vector_id, :similar_gi_vector_id], :name => "gi1_gi2"
    add_index :gi_vector_similarities, [:similar_gi_vector_id, :gi_vector_id], :name => "gi2_gi1"
  end

  def self.down
    drop_table :gi_vector_similarities
  end
end
