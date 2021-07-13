from flask import Flask, render_template, session
from flask_bootstrap import Bootstrap
from flask_wtf import Form
from wtforms import StringField, SubmitField
from wtforms.validators import Required, Length
app = Flask(__name__)
app.config['SECRET_KEY'] = 'top secret!'
bootstrap = Bootstrap(app)

class NameForm(Form):
    name = StringField('What is your name?', validators=[Required(),
                                                         Length(1, 16)])
    submit = SubmitField('Submit')


@app.route('/', methods=['GET', 'POST'])
def index():
    print("IN INDEX1")
    name = None
    form = NameForm()
    if form.validate_on_submit():
        name = form.name.data
        form.name.data = ''
        print(name)
        if name.startswith('n'):
            session['count'] = [1, 2, 3]
        else:
            session['count'] = [4, 5, 6]
    print("IN INDEX2",name,session['count'])
    return render_template('index.html', form=form, name=name,count=session['count'])

@app.errorhandler(500)
def internal_server_error(error):
    app.logger.error('Server Error: %s', (error))
    return render_template('500.htm'), 500

if __name__ == '__main__':
   app.run()
   app.config['PROPAGATE_EXCEPTIONS'] = True