from flask import Flask, redirect, url_for, request, session, render_template
from autoCurator import curate, nlp
import time
app = Flask(__name__)


@app.route('/')
def index():
    return render_template('home.html')

@app.route('/home', methods=['POST', 'GET'])
def set_session():
    session.clear()
    if request.method == 'POST':
        user = request.form['nm']
        session['name'] = request.form['nm']
        session['name'] = session.get('name')
        session['year1'] = request.form['year1']
        session['year1'] = session.get('year1')
        #print(session['year1'])

        session['year2'] = request.form['year2']
        session['year2'] = session.get('year2')
        #print(session['year2'])

        session['NumOfPub'] = request.form['NumOfPub']
        session['NumOfPub'] = session.get('NumOfPub')

        #print(session['NumOfPub'])
        return redirect(url_for('results'))
    return render_template("home.html")

@app.route('/resultss', methods=['POST', 'GET'])
def resultss():
    session['name']= request.args.get('name', None)
    cPapers,nAbs,nSP = curate(session.get('name'), session['year1'], session['year2'], session['NumOfPub'])
    #print(len(cPapers))
    results = formulateSNP(cPapers)
    session['results'] = results[0]
    session['PaperCount'] = results[1]
    session['SNPcount'] = results[2]
    session['AbstractCnt'] = nAbs
    session['AbstractSNPCnt'] = nSP
    print("HERE", session['AbstractCnt'], session['AbstractSNPCnt'])
    print("results ", session.get('SNPcount')," from ",session.get('PaperCount'), "Papers" )
    return render_template("results.html",disease=session.get('name'), rows=session.get('results'),
                           SNPCount=session.get('SNPcount'), PaperCount=session.get('PaperCount'),AbsrtactCount=session.get('AbstractCnt'),AbsrtactSNPCount=session.get('AbstractSNPCnt'))

@app.route('/results', methods=['POST', 'GET'])
def results():
    cPapers, nAbs, nSP= curate(session.get('name'), session['year1'], session['year2'], session['NumOfPub'])
    ##print(len(cPapers))
    results = formulateSNP(cPapers)
    session['results'] = results[0]
    print(results[0])
    session['PaperCount'] = results[1]
    session['SNPcount'] = results[2]
    session['AbstractCnt']=nAbs
    session['AbstractSNPCnt'] = nSP
    print("HERE",session['AbstractCnt'],session['AbstractSNPCnt'])
    print("results ", session.get('SNPcount')," from ",session.get('PaperCount'), "Papers" )
    return render_template("results.html",disease=session.get('name'), rows=session.get('results'),
                           SNPCount=session.get('SNPcount'), PaperCount=session.get('PaperCount'),AbsrtactCount=session.get('AbstractCnt'),AbsrtactSNPCount=session.get('AbstractSNPCnt'))

@app.route('/About')
def GoToAbout():
    return render_template('about.html')



def formulateSNP(d):
    data = []
    countPapers = 0
    for i in d:
        if i[8]:
            countPapers = countPapers + 1
        for n in (n for n in i[8] if (n[1] != '' or n[2] != '')):
            #print(i[0], n[0], n[1])
            try:
                x = float(n[1])
                n[1] = x
            except:
                pass
            try:
                x = int(i[5])
                i[5] = x
            except:
                pass
            try:
                x = int(i[6])
                i[6] = x
            except:
                pass
            finally:
                data.append([i[0], i[1], i[2], n[0], n[1], n[2], n[3], i[4], i[5], i[6],i[9]])
    AllSnp=[]
    for i in data:
        AllSnp.append(i[3])
        #print(i)

    for i in data:
        #print(i[3])
        #print(AllSnp.count(i[3]))
        i[10]=AllSnp.count(i[3])

    keys = ["PMID", "Title", "Date", "SNP", "Pvalue", "ORvalue", "EvidenceSent", "Ethnicity", "PatientSize","ControlSize","Frequency"]
    countSNP = len(data)

    Results = dict(zip(keys, zip(*data)))
    #print(Results)
    R = [Results, countPapers, countSNP]
    return R



@app.errorhandler(500)
def internal_server_error(error):
    app.logger.error('Server Error: %s', (error))
    return render_template('500.htm'), 500


app.secret_key = 'tsdhisiusdfdsfaSecsdfsdfrfghdetkey'
if __name__ == "__main__":
    st = time.time()
    app.run(threaded=True, debug=True)
    f.write("Run time in main",time.time()  - st)
