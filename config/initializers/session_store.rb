# Be sure to restart your server when you modify this file.

# Your secret key for verifying cookie session data integrity.
# If you change this key, all old sessions will become invalid!
# Make sure the secret is at least 30 characters and all random, 
# no regular words or you'll be exposed to dictionary attacks.
ActionController::Base.session = {
  :key         => '_mink_session',
  :secret      => '96510e8531bb1858a069fb3a4fabfd6f4077faba668d7a506edc56dc01e0a17eb93e3819788668942cffa7a17b75ac2450bced97c22dddad772b7c6ce3479b37'
}

# Use the database for sessions instead of the cookie-based default,
# which shouldn't be used to store highly confidential information
# (create the session table with "rake db:sessions:create")
# ActionController::Base.session_store = :active_record_store
