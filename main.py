from flask import Flask, request
from flask_pymongo import PyMongo
from bson.objectid import ObjectId
import smtplib
import ssl

port = 465
context = ssl.create_default_context()

app = Flask(__name__)
app.config["MONGO_URI"] = "mongodb://localhost:27017/sigtime"
mongo = PyMongo(app)

db_operations = mongo.db.results

@app.route('/v1/analysis', methods=['POST'])
def analysis():
    req_data = request.get_json(force=True)
    id = req_data['id']

    # todo 여기다가 작업코드 넣어

    receiver_email = db_operations.find_one({ "_id": ObjectId(id)})["email"]
    print(receiver_email)
    result = db_operations.update_one({ "_id": ObjectId(id) }, {"$set": { "complete": True }})
    if result is not None:
        with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
            to = [receiver_email]
            subject = 'test url'
            email_text = """\
From: %s
To: %s
Subject: %s

http://localhost:3000/result/%s
            """ % ("q1w2e30630@gmail.com", ", ".join(to), subject, id)
            print(email_text)
            server.login("q1w2e30630@gmail.com", "cjswo0630135")
            server.sendmail("q1w2e30630@gmail.com", receiver_email, email_text)
    return 'complete'

app.run(port=5001)