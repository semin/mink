# Add your own tasks in files placed in lib/tasks ending in .rake,
# for example lib/tasks/capistrano.rake, and they will automatically be available to Rake.

require(File.join(File.dirname(__FILE__), 'config', 'boot'))

require 'rake'
require 'rake/testtask'
require 'rake/rdoctask'

require 'tasks/rails'
require 'thinking_sphinx/tasks'

include FileUtils

#RakeFileUtils.verbose(false)

Rake.application.options.trace = true

#Rake.application.instance_variable_get(:@tasks).delete("db:schema:dump")
#namespace(:db) { namespace(:schema) { task(:dump) { puts "Schema dump disabled" } } }

class RakeLogger < Logger
  def format_message(severity, timestamp, progname, msg)
          "[#{timestamp.to_formatted_s(:db)} #{severity}] #{msg}\n"
  end
end

$logger       = RakeLogger.new(STDOUT)
$logger.level = RakeLogger::DEBUG
