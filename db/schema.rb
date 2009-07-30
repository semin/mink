# This file is auto-generated from the current state of the database. Instead of editing this file, 
# please use the migrations feature of Active Record to incrementally modify your database, and
# then regenerate this schema definition.
#
# Note that this schema.rb definition is the authoritative source for your database schema. If you need
# to create the application database on another system, you should be using db:schema:load, not running
# all the migrations from scratch. The latter is a flawed and unsustainable approach (the more migrations
# you'll amass, the slower it'll run and the greater likelihood for issues).
#
# It's strongly recommended to check this file into your version control system.

ActiveRecord::Schema.define(:version => 0) do

  # 'scop' table
  create_table :scop, :force => true do |t|
    t.belongs_to  :parent
    t.integer     :lft
    t.integer     :rgt
    t.string      :type
    t.integer     :sunid
    t.string      :stype
    t.string      :sccs
    t.string      :sid
    t.string      :description
  end

  add_index :scop, :sunid
  add_index :scop, :parent_id
  add_index :scop, :lft
  add_index :scop, :rgt
  add_index :scop, [:id, :type]


  # 'mink_vectors' table
  create_table :mink_vectors, :force => true do |t|
    t.belongs_to  :scop
    t.string      :sid
    t.integer     :sunid
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
  end

  add_index :mink_vectors, :sid,    :unique => true
  add_index :mink_vectors, :sunid,  :unique => true


  # 'norm_mink_vectors' table
  create_table :norm_mink_vectors, :force => true do |t|
    t.belongs_to  :scop
    t.belongs_to  :mink_vector
    t.string      :sid
    t.integer     :sunid
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
  end

  add_index :norm_mink_vectors, :sid,    :unique => true
  add_index :norm_mink_vectors, :sunid,  :unique => true

end
