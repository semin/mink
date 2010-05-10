class CreateMinkSearchHits < ActiveRecord::Migration
  def self.up
    create_table :mink_search_hits, :force => true do |t|
      t.belongs_to  :mink_search
      t.belongs_to  :norm_mink_vector
      t.float       :distance
    end

    add_index :mink_search_hits, :distance
  end

  def self.down
    drop_table :mink_search_hits
  end
end
