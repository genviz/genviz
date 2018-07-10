from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth import get_user_model
User = get_user_model()
from app.models import Patient


class SignUpForm(UserCreationForm):
    first_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    last_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    email = forms.EmailField(max_length=254, help_text='Required. Inform a valid email address.')

    class Meta:
        model = User
        fields = ('first_name', 'last_name', 'email', 'password1', 'password2', )

class PatientForm(forms.ModelForm):
    class Meta:
        model = Patient
        fields = (
        	'identifier',
        	'first_name',
        	'last_name',
        	'sex',
        	'birthday',
        	'sample_date',
        	'phone',
        	'diagnosis',
        	'email'
        )
        widgets = {
            'birthday': forms.DateInput(attrs={'class':'datepicker'}),
            'sample_date': forms.DateInput(attrs={'class':'datepicker'}),
        }

    def clean(self):
        cleaned_data = super().clean()
        birthday = cleaned_data.get("birthday")
        sample_date = cleaned_data.get("sample_date")
        if sample_date < birthday:
            msg = u"Sample date should be greater than Birthday."
            self._errors["sample_date"] = self.error_class([msg])
        return cleaned_data
