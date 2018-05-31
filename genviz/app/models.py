
from django.contrib.auth.models import AbstractUser, BaseUserManager
from django.db import models
from django.utils.translation import ugettext_lazy as _


class UserManager(BaseUserManager):
    use_in_migrations = True

    def _create_user(self, email, password, **extra_fields):
        if not email:
            raise ValueError('The given email must be set')
        email = self.normalize_email(email)
        user = self.model(email=email, **extra_fields)
        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_user(self, email, password=None, **extra_fields):
        extra_fields.setdefault('is_staff', False)
        extra_fields.setdefault('is_superuser', False)
        return self._create_user(email, password, **extra_fields)

    def create_superuser(self, email, password, **extra_fields):
        extra_fields.setdefault('is_staff', True)
        extra_fields.setdefault('is_superuser', True)

        if extra_fields.get('is_staff') is not True:
            raise ValueError('Superuser must have is_staff=True.')
        if extra_fields.get('is_superuser') is not True:
            raise ValueError('Superuser must have is_superuser=True.')

        return self._create_user(email, password, **extra_fields)


class Patient(models.Model):
    first_name = models.CharField(max_length=100)
    last_name  = models.CharField(max_length=100)
    identifier = models.CharField(max_length=30, unique=True)

    def __str__(self):
        return self.full_name()

    def full_name(self):
        return "%s %s" % (self.first_name, self.last_name)


class User(AbstractUser):
    username     = None
    email        = models.EmailField(_('email address'), unique=True)
    is_doctor    = models.BooleanField('student status', default=False)
    is_biologist = models.BooleanField('teacher status', default=True)
    patients     = models.ManyToManyField(Patient)

    USERNAME_FIELD = 'email'
    REQUIRED_FIELDS = []

    objects = UserManager()

class Annotation(models.Model):
    OPERATION_CHOICES = (
        ('ins', 'Insertion'),
        ('del', 'Deletion'),
        ('insdel', 'Insertion-Deletion'),
        ('sub', 'Substitution'),
        ('dup', 'Duplication'),
    )
    author    = models.ForeignKey(User, on_delete=models.CASCADE)
    seq_id    = models.CharField(max_length=30)
    after     = models.IntegerField(null=True)
    start     = models.IntegerField(null=True)
    end       = models.IntegerField(null=True)
    sequence  = models.CharField(max_length=1000)
    operation = models.CharField(max_length=10, choices=OPERATION_CHOICES)
    comment   = models.TextField(null=True)
    patient   = models.ForeignKey(Patient, on_delete=models.CASCADE)
    source    = models.CharField(max_length=30, default='manual')