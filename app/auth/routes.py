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

        user = select(u for u in db.User if form.username.data == u.username).first()
        if user is None or not user.check_password(form.password.data):
            flash('Invalid Username or Password')
            return redirect(url_for('auth.auth'))
        login_user(user)
        next_page = request.args.get('next')
        if not next_page or url_parse(next_page).netloc != '':

            next_page = url_for('main.home')
        return redirect(next_page)
    return render_template('index.html', form=form)


@auth_bp.route('/auth', methods=['GET', 'POST'])
def auth():
    """
    Handles user authentication process by logging out the current user and
    redirecting to the authentication index.

    :raises RuntimeError: If called outside an application context.
    :return: A redirect response to the authentication index page.
    :rtype: Response
    """
    logout_user()
    return redirect(url_for('auth.index'))


@auth_bp.route('/signin', methods=['GET', 'POST'])
def signin():
    """
    Renders the sign-in page for the application.

    This function handles the HTTP GET and POST methods for the sign-in page. It
    is typically used to either display the sign-in form to the user or handle the
    form submission for user authentication.

    :return: A rendered HTML page for the sign-in form.
    :rtype: str
    """
    return render_template('signin.html')


@auth_bp.route('/register', methods=['GET', 'POST'])
def register():  # the view function that handles user registrations
    """
    Handles user registration by rendering a registration form and processing user
    input from the form. If the form passes validation, a new user is created and
    saved to the database. Upon successful registration, the user is redirected
    to the login page, and a success message is displayed.

    :returns: If the form is submitted and validated, redirects the user to the
        login page. Otherwise, renders the registration page with the form.
    :rtype: werkzeug.wrappers.Response or flask.Response
    """
    form = RegistrationForm()  # instantiates an object of RegistrationForm
    if form.validate_on_submit():
        # Creates a user and commits to the database
        db.User(username=form.username.data, email=form.email.data, fullname=form.fullname.data,
                password_hash=db.User.set_password(form.password.data))
        flash('Congratulations, you are now a registered user!')  # flashes the success message
        return redirect(url_for('auth.auth'))  # redirects to the login page
    return render_template('register.html', form=form)  # renders the registration template
