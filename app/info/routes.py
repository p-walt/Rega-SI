from info import info_bp
from flask import render_template


@info_bp.route('/contact', methods=['GET', 'POST'])
def contact():
    return render_template('contact.html')


@info_bp.route('/about', methods=['GET', 'POST'])
def about():
    return render_template('about.html')


@info_bp.route('/theory', methods=['GET', 'POST'])
def theory():
    return render_template('theory.html')


@info_bp.route('/references', methods=['GET', 'POST'])
def references():
    return render_template('references.html')


@info_bp.route('/user_guide', methods=['GET', 'POST'])
def user_guide():
    return render_template('user_guide.html')


