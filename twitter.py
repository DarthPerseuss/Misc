# This program prints tweets that contain a certain piece of string
# set by 'track_string' variable

import tweepy
from tweepy.streaming import StreamListener
from tweepy import OAuthHandler
from tweepy import Stream
import json
import numpy as np

# Variables that contains the user credentials to access Twitter API
access_token = "909412155144814593-zUs2ce9AhWLpFQEmJqyl2fGoJ72ulss"
access_token_secret = "NbZa9fwNS5Y1lkhqxRc62yjwtFx5kCcNowXUNRgV82BEi"
consumer_key = "xqJdht7Y1OQtVDCKxj3BUuDYk"
consumer_secret = "TvV3AqZrAViJKunICtWoymKsvpQenDqyi1XX7c0TyFiSo3NPd1"

# This handles Twitter authentication and the connection to Twitter Streaming API

auth = OAuthHandler(consumer_key, consumer_secret)
auth.set_access_token(access_token, access_token_secret)

api = tweepy.API(auth)

# Print the tweets in my homepage
# public_tweets = api.home_timeline()
# for tweet in public_tweets:
#     print(tweet.text)


# This listener will print out all Tweets it receives
class PrintListener(tweepy.StreamListener):
    def on_data(self, data):
        # Decode the JSON data
        tweet = json.loads(data)

        # Print out the Tweet
        print('@%s: %s' % (tweet['user']['screen_name'], tweet['text'].encode('ascii', 'ignore')))

    def on_error(self, status):
        print(status)

    
if __name__ == '__main__':
    listener = PrintListener()
    track_string = 'Islamic Republic'
    # Show system message
    print('I will now print Tweets containing {}! ==>'.format(track_string))

    # Authenticate
    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
    auth.set_access_token(access_token, access_token_secret)

    # Connect the stream to our listener
    stream = tweepy.Stream(auth, listener)
    stream.filter(track=[track_string])

