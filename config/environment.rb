# Be sure to restart your server when you modify this file

# Specifies gem version of Rails to use when vendor/rails is not present
RAILS_GEM_VERSION = '2.3.8' unless defined? RAILS_GEM_VERSION

# Bootstrap the Rails environment, frameworks, and default configuration
require File.join(File.dirname(__FILE__), 'boot')

Rails::Initializer.run do |config|
  # Settings in config/environments/* take precedence over those specified here.
  # Application configuration should go into files in config/initializers
  # -- all .rb files in that directory are automatically loaded.

  # Add additional load paths for your own custom dirs
  # config.load_paths += %W( #{RAILS_ROOT}/extras )

  # Specify gems that this application depends on and have them installed with rake gems:install
  # config.gem "bj"
  # config.gem "hpricot", :version => '0.6', :source => "http://code.whytheluckystiff.net"
  # config.gem "sqlite3-ruby", :lib => "sqlite3"
  # config.gem "aws-s3", :lib => "aws/s3"

  config.gem 'uuid'
  config.gem 'narray'
  config.gem 'paperclip'
  config.gem 'configatron'
  config.gem 'daemon-spawn', :lib => 'daemon_spawn'
  config.gem 'ar-extensions'
  config.gem 'parallel-forkmanager', :lib => 'forkmanager'
  config.gem 'thinking-sphinx', :lib => 'thinking_sphinx', :version => '1.3.16'

  # Only load the plugins named here, in the order given (default is alphabetical).
  # :all can be used as a placeholder for all plugins not explicitly named
  # config.plugins = [ :exception_notification, :ssl_requirement, :all ]

  # Skip frameworks you're not going to use. To use Rails without a database,
  # you must remove the Active Record framework.
  # config.frameworks -= [ :active_record, :active_resource, :action_mailer ]

  # Activate observers that should always be running
  # config.active_record.observers = :cacher, :garbage_collector, :forum_observer

  # Set Time.zone default to the specified zone and make Active Record auto-convert to this zone.
  # Run "rake -D time" for a list of tasks for finding time zone names.
  config.time_zone = 'UTC'

  # The default locale is :en and all translations from config/locales/*.rb,yml are auto loaded.
  # config.i18n.load_path += Dir[Rails.root.join('my', 'locales', '*.{rb,yml}')]
  # config.i18n.default_locale = :de

  config.after_initialize do
    configatron.username            = '26e60a40d5d5e3d64afd5058600bc67f9e2714c0'
    configatron.password            = '7dad210b9683a031ab75ce48ce3f40dafbc44c3a'

    configatron.scop_dir            = Pathname.new("/merlin/Store/SCOP")

    configatron.bin_dir             = Rails.root.join("bin")
    configatron.mink                = configatron.bin_dir.join("mink")
    configatron.minkproc            = configatron.bin_dir.join("minkproc")
    configatron.batch_mink          = configatron.bin_dir.join("batch_mink.py")
    configatron.proc_mink           = configatron.bin_dir.join("proc_mink.py")
    configatron.gi                  = configatron.bin_dir.join("GI")

    configatron.pub_dir             = Rails.root.join("public")
    configatron.mink_dir            = configatron.pub_dir.join("mink")
    configatron.fig_dir             = configatron.pub_dir.join("figures")
    configatron.scop_fig_dir        = configatron.fig_dir.join("scop")
    configatron.mink_scop_dir       = configatron.pub_dir.join("scop")
    configatron.mink_scop_gi_dir    = configatron.mink_scop_dir.join("gi")
    configatron.mink_scop_git_dir   = configatron.mink_scop_dir.join("git")
    configatron.mink_scop_pdb_dir   = configatron.mink_scop_dir.join("pdb")
    configatron.mink_scop_mink_dir  = configatron.mink_scop_dir.join("mink")

    configatron.src_dir                 = Rails.root.join("src")
    configatron.bin_dir                 = Rails.root.join("bin")
    configatron.calculate_distances_src = configatron.src_dir.join("calculate_distances.c")
    configatron.calculate_distances_bin = configatron.bin_dir.join("calculate_distances")
    configatron.gi_src                  = configatron.src_dir.join("GI.c")
    configatron.gi_bin                  = configatron.bin_dir.join("GI")
    configatron.git_src                 = configatron.src_dir.join("GIT.c")
    configatron.git_bin                 = configatron.bin_dir.join("GIT")
    configatron.avg_gauss_tbl_smt_rep   = configatron.mink_scop_dir.join("AverageGaussTableSmoothRepresentation")

    configatron.tmp_dir               = Rails.root.join("tmp")
    configatron.mink_vectors_csv      = configatron.tmp_dir.join("mink_vectors.csv")
    configatron.norm_mink_vectors_csv = configatron.tmp_dir.join("norm_mink_vectors.csv")
    configatron.gi_vectors_csv        = configatron.tmp_dir.join("gi_vectors.csv")
    configatron.git_vectors_csv       = configatron.tmp_dir.join("git_vectors.csv")
    configatron.norm_gi_vectors_csv   = configatron.tmp_dir.join("norm_gi_vectors.csv")
    configatron.norm_git_vectors_csv  = configatron.tmp_dir.join("norm_git_vectors.csv")
    configatron.mink_vector_similarities_csv      = configatron.tmp_dir.join("mink_vector_similarities.csv")
    configatron.norm_mink_vector_similarities_csv = configatron.tmp_dir.join("norm_mink_vector_similarities.csv")
    configatron.gi_vector_similarities_csv        = configatron.tmp_dir.join("gi_vector_similarities.csv")
    configatron.git_vector_similarities_csv       = configatron.tmp_dir.join("git_vector_similarities.csv")
    configatron.norm_gi_vector_similarities_csv   = configatron.tmp_dir.join("norm_gi_vector_similarities.csv")
    configatron.norm_git_vector_similarities_csv  = configatron.tmp_dir.join("norm_git_vector_similarities.csv")

    require 'mink/distance'
    require 'gi/distance'
    require 'mink_search_job'
    require_dependency "scop"

  end
end
