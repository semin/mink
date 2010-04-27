class GiSearch < ActiveRecord::Base

  include Gi::Distance
  include FileUtils

  has_many :gi_search_hits

  has_attached_file :pdb

  validates_attachment_presence :pdb,
                                :message => "file must be provided."

  validates_attachment_size :pdb,
                            :in => 1..1.megabyte,
                            :message => "must be in range of 1 byte to 1 megabyte."

  validates_numericality_of :cutoff,
                            :greater_than_or_equal_to => 0,
                            :less_than_or_equal_to => 1

  def to_param
    self.uuid
  end

  def generate_images
    cwd     = pwd
    pdbpath = Pathname.new(pdb.path)

    chdir pdbpath.dirname

    pdb       = pdbpath.basename
    stem      = pdbpath.basename(pdbpath.extname)
    input     = "#{stem}.molinput"
    fig500    = "#{stem}-500.png"
    fig100    = "#{stem}-100.png"
    mol       = `molauto -notitle -nice #{pdb}`.split("\n")
    mol[5,0]  = "  background grey 1;"

    File.open(input, "w") { |f| f.puts mol.join("\n") }
    system "molscript -r < #{input} | render -png #{fig500} -size500x500"
    system "convert #{fig500} -resize 100x100 #{fig100}"
    chdir cwd
  end

  def big_image
    pdbpath = Pathname.new(pdb.path)
    "/" + pdbpath.dirname.join(pdbpath.basename(pdbpath.extname).to_s + "-500.png").relative_path_from(Rails.root.join("public")).to_s
  end

  def small_image
    pdbpath = Pathname.new(pdb.path)
    "/" + pdbpath.dirname.join(pdbpath.basename(pdbpath.extname).to_s + "-100.png").relative_path_from(Rails.root.join("public")).to_s
  end

end
