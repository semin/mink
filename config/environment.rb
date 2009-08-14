# Be sure to restart your server when you modify this file

# Specifies gem version of Rails to use when vendor/rails is not present
RAILS_GEM_VERSION = '2.3.3' unless defined? RAILS_GEM_VERSION

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

  config.gem 'narray'
  config.gem 'uuid'
  config.gem 'configatron'
  config.gem 'freelancing-god-thinking-sphinx', :lib => 'thinking_sphinx', :source => 'http://gems.github.com'
  config.gem 'alexvollmer-daemon-spawn', :lib => 'daemon-spawn', :source => 'http://gems.github.com'
  config.gem 'thoughtbot-paperclip', :lib => 'paperclip', :source => 'http://gems.github.com'

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

  # config ActionMailer
  config.action_mailer.delivery_method = :sendmail
  config.action_mailer.smtp_settings = {
    :address        => "localhost",
    :port           => 25,
    :domain         => "cryst.bioc.cam.ac.uk",
    :authentication => :plain,
  }

  config.after_initialize do
    configatron.username      = '26e60a40d5d5e3d64afd5058600bc67f9e2714c0'
    configatron.password      = '7dad210b9683a031ab75ce48ce3f40dafbc44c3a'
    configatron.scop_pdb_dir  = Pathname.new("/BiO/Store/SCOP/pdbstyle")
    configatron.mink          = Rails.root.join("bin", "mink")
    configatron.minkproc      = Rails.root.join("bin", "minkproc")
    configatron.mink_dir      = Rails.root.join("public", "minkscop")
    configatron.fig_dir       = Rails.root.join("public", "figures")
    configatron.scop_fig_dir  = configatron.fig_dir.join("scop")

    require 'mink/distance'
    require 'mink_search_job'
    require_dependency "scop"
  end
end
