from wtforms import StringField, SubmitField, SelectField, PasswordField, BooleanField, validators
from wtforms.validators import DataRequired, Email, EqualTo, ValidationError
from flask_wtf import FlaskForm
from flask_wtf.file import FileField
from pony.orm import select
from app import db


class SubmitForm(FlaskForm):  # this class defines the login form fields
    name = StringField('Name of Compound')
    smiles = StringField('SMILES string of Compound')
    radical = SelectField('Radical to add (cf3, cf2h, ipr)', choices=["-select-", "cf3", "cf2h", "ipr"])
    hpc_calcs = BooleanField('Calculate accurate energies using the HPC?')
    file = FileField('File')
    submit_single = SubmitField(render_kw={"onclick": "loading()", "class": "btn btn-primary"})
    submit_multiple = SubmitField(render_kw={"onclick": "loading()", "class": "btn btn-primary"})


class LoginForm(FlaskForm):
    username = StringField('Username', [DataRequired()], render_kw={"placeholder": "Username",
                                                                    "class": "form-control form-control-lg"})
    password = PasswordField('Password', [DataRequired()], render_kw={"placeholder": "Password",
                                                                      "class": "form-control form-control-lg"})
    submit = SubmitField(render_kw={"class": "btn btn-primary"})


class RegistrationForm(FlaskForm):  # this class defines the registration form fields
    username = StringField('Username', validators=[DataRequired()], render_kw={"placeholder": "Username",
                                                                               "class": "form-control form-control-lg"})
    email = StringField('Email', validators=[DataRequired(), Email()], render_kw={"placeholder": "Email",
                                                                                  "class": "form-control form-control-lg"})
    fullname = StringField('Full Name', validators=[DataRequired()], render_kw={"placeholder": "Full Name",
                                                                                "class": "form-control form-control-lg"})
    password = PasswordField('Password', validators=[DataRequired()], render_kw={"placeholder": "Password",
                                                                                 "class": "form-control form-control-lg"})
    password2 = PasswordField(
        'Repeat Password', validators=[DataRequired(), EqualTo('password')], render_kw={"placeholder": "Repeat Password",
                                                                                        "class": "form-control form-control-lg"})
    # a user is asked to type the password twice to reduce the risk of a typo
    submit = SubmitField('Register', render_kw={"class": "btn btn-primary"})

    """WTForms takes methods like validate_<field_name> as custom validators 
    and invokes them in addition to the stock validators. The validators
    validate_username and validate_email are used to make sure that the 
    username and email address entered by the user are not already in 
    the database, so these two methods issue database queries expecting 
    there will be no results. In the event a result exists, a validation 
    error is triggered by raising ValidationError. The message included 
    as the argument in the exception will be the message that will be 
    displayed next to the field for the user to see."""
    """The database queries are now using select() from Pony.orm""" # TODO RA: Should be a comment. 
    def validate_username(self, username):
        user = select(u for u in db.User if username.data == u.username).first()
        if user is not None:
            raise ValidationError('This username is already taken.')

    def validate_email(self, email):
        user = select(u for u in db.User if email.data == u.email).first()
        if user is not None:
            raise ValidationError('There is already an account associated with this email.')
