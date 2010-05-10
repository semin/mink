class CreateNews < ActiveRecord::Migration
  def self.up
    create_table :news, :force => true do |t|
      t.date    :date
      t.string  :title
      t.text    :content
    end

    add_index :news, :date
  end

  def self.down
    drop_table :news
  end
end
