class CreateScops < ActiveRecord::Migration
  def self.up
    create_table :scops, :force => true do |t|
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

    add_index :scops, :sunid
    add_index :scops, :parent_id
    add_index :scops, :lft
    add_index :scops, :rgt
    add_index :scops, [:id, :type]
  end

  def self.down
    drop_table :scops
  end
end
