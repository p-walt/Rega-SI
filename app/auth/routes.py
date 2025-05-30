from auth import auth_bp
from flask_login import login_user, logout_user
from flask import render_template, request, url_for, redirect, flash
from werkzeug.urls import url_parse
from forms import LoginForm, RegistrationForm
from pony.orm import select
from app import db


# initial log in page
@auth_bp.route('/login', methods=['GET', 'POST'])
def index():
    form = LoginForm()
    if request.method == 'POST':
        '''The form.validate_on_submit returns True when the browser sends the POST 
                request as a result of the user pressing the submit button and if all the fields
                passes validation. It returns False when the browser sends the GET request to 
                receive the web page with the form or if at least one field fails validation.'''
        user = select(u for u in db.User if form.username.data == u.username).first()
        '''The select function will search through all of the User entities in the
        database and will return a query that only includes the objects that have 
        a matching username. Since there is only one or zero results, the query is 
        completed by calling first(), which will return the user object if it exists, 
        or None if it does not.'''
        if user is None or not user.check_password(form.password.data):
            '''If it got a match for the username that was provided, it can next check 
            if the password came with the form is valid. This is done by invoking the 
            check_password() method defined in models.py. This will take the password 
            hash stored with the user and determine if the password entered in the 
            form matches the hash or not. In either of two possible error conditions -
            the invalid username or the incorrect password - the error message is
            flashed, and the user is redirected back to the login prompt to try again.'''
            flash('Invalid Username or Password')
            return redirect(url_for('auth.auth'))
        login_user(user)
        '''If the username and password are both correct, then the login_user() function
        from Flask-Login is called. This function will register the user as logged in, 
        which means that any future pages the user navigates to will have the current_user 
        variable set to that user.'''
        next_page = request.args.get('next')
        '''The next query string argument is set to the original URL, 
        so the application can use that to redirect back after login.'''
        if not next_page or url_parse(next_page).netloc != '':
            '''If the login URL does not have a next argument or the 
            next argument is set to a full URL that includes a domain 
            name, then the user is redirected to the index page.'''
            next_page = url_for('main.home')
            '''If the login URL includes a next argument that is set 
            to a relative path (a URL without the domain portion), 
            then the user is redirected to that URL.'''
        return redirect(next_page)
    return render_template('index.html', form=form)


@auth_bp.route('/auth', methods=['GET', 'POST'])
def auth():
    logout_user()
    return redirect(url_for('auth.index'))


@auth_bp.route('/signin', methods=['GET', 'POST'])
def signin():
    return render_template('signin.html')


@auth_bp.route('/register', methods=['GET', 'POST'])
def register():  # the view function that handles user registrations
    form = RegistrationForm()  # instantiates an object of RegistrationForm
    if form.validate_on_submit():
        '''The form.validate_on_submit returns True when the browser sends the POST
        request as a result of the user pressing the submit button and if all the fields
        passes validation. It returns False when the browser sends the GET request to
        receive the web page with the form or if at least one field fails validation.'''
        # Creates a user and commits to the database
        db.User(username=form.username.data, email=form.email.data, fullname=form.fullname.data,
                password_hash=db.User.set_password(form.password.data))
        flash('Congratulations, you are now a registered user!')  # flashes the success message
        return redirect(url_for('auth.auth'))  # redirects to the login page
    return render_template('register.html', form=form)  # renders the registration template
