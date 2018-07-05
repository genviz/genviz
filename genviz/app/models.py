from django.contrib.auth.models import AbstractUser, BaseUserManager
from django.db import models
from django.utils.translation import ugettext_lazy as _
import hgvs.location
import hgvs.posedit

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
    SEX_CHOICES = (
        ('F','Female'),
        ('M','Male')
    )

    first_name = models.CharField(max_length=100) 
    last_name  = models.CharField(max_length=100)
    identifier = models.CharField(max_length=30, unique=True)
    age = models.IntegerField()
    sex = models.CharField(max_length=1, choices=SEX_CHOICES, default='M')
    phone  = models.CharField(max_length=25, default=None,null=True)
    diagnosis = models.CharField(max_length=100) 
    email = models.EmailField(max_length=100)
   
    def __str__(self):
        return self.full_name()

    def full_name(self):
        return "%s %s" % (self.first_name, self.last_name)




class User(AbstractUser):
    username     = models.CharField(max_length=100, null=True)
    email        = models.EmailField(_('email address'), unique=True)
    is_doctor    = models.BooleanField('student status', default=False)
    is_biologist = models.BooleanField('teacher status', default=True)
    patients = models.ManyToManyField(Patient)

    USERNAME_FIELD = 'email'
    REQUIRED_FIELDS = []

    objects = UserManager()

    def full_name(self):
        return self.first_name + ' ' + self.last_name



class Variation(models.Model):
    OPERATION_CHOICES = (
        ('ins', 'Insertion'),
        ('del', 'Deletion'),
        ('insdel', 'Insertion-Deletion'),
        ('delins', 'Deletion-Insertion'),
        ('sub', 'Substitution'),
        ('dup', 'Duplication'),
    )
    author    = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    acc_id    = models.CharField(max_length=30)
    start     = models.IntegerField(null=True)
    end       = models.IntegerField(null=True)
    ref       = models.CharField(max_length=1000)
    alt       = models.CharField(max_length=1000)
    operation = models.CharField(max_length=10, choices=OPERATION_CHOICES)
    comment   = models.TextField(null=True)
    patient   = models.ForeignKey(Patient, on_delete=models.CASCADE, null=True)
    source    = models.CharField(max_length=30, default='Unknown')
    url       = models.CharField(max_length=255, null=True)
    coordinate_type = models.CharField(max_length=1, default='n')

    def hgvs(self):
        if self.operation == 'ins':
            return '{acc}:{coordinate_type}.{start}_{end}ins{alt}'.format(
                acc=self.acc_id,
                coordinate_type=self.coordinate_type,
                start=self.start,
                end=self.end,
                alt=self.alt
            )
        elif self.operation == 'sub':
            return '{acc}:{coordinate_type}.{loc}{ref}>{alt}'.format(
                acc=self.acc_id,
                coordinate_type=self.coordinate_type,
                loc=self.start,
                ref=self.ref,
                alt=self.alt
            )
        elif self.operation == 'delins':
            if self.start - self.end + 1 > 1:
                return '{acc}:{coordinate_type}.{start}_{end}delins{alt}'.format(
                    acc=self.acc_id,
                    coordinate_type=self.coordinate_type,
                    start=self.start,
                    end=self.end,
                    alt=self.alt
                )
            else:
                return '{acc}:{coordinate_type}.{start}delins{alt}'.format(
                    acc=self.acc_id,
                    coordinate_type=self.coordinate_type,
                    start=self.start,
                    alt=self.alt
                )
        elif self.operation == 'del':
            if self.start - self.end + 1 > 1:
                return '{acc}:{coordinate_type}.{start}_{end}del'.format(
                    acc=self.acc_id,
                    coordinate_type=self.coordinate_type,
                    start=self.start,
                    end=self.end
                )
            else:
                return '{acc}:{coordinate_type}.{start}del'.format(
                    acc=self.acc_id,
                    coordinate_type=self.coordinate_type,
                    start=self.start
                )
        elif self.operation == 'dup':
            if self.start - self.end + 1 > 1:
                return '{acc}:{coordinate_type}.{start}_{end}dup'.format(
                    acc=self.acc_id,
                    coordinate_type=self.coordinate_type,
                    start=self.start,
                    end=self.end
                )
            else:
                return '{acc}:{coordinate_type}.{start}dup'.format(
                    acc=self.acc_id,
                    coordinate_type=self.coordinate_type,
                    start=self.start
                )

    @staticmethod
    def from_hgvs(variation, source):
        """
        Parses a string that describes a change in hgvs format and returns Variation object
        """
        hgvs_variation = hgvs.parse_hgvs_variant(variation)
        return self.from_hgvs_obj(hgvs_variation, source)

    @staticmethod
    def from_hgvs_obj(hgvs_variation, source):
        """
        Maps a hgvs.SequenceVariant to a Variation object
        """
        pos = hgvs_variation.posedit.pos
        start, end = pos.start, pos.end
        acc_id = hgvs_variation.ac
        edit = hgvs_variation.posedit.edit
        ref, alt = edit.ref, edit.alt

        change = str(edit)
        if '>' in change:
            operation = 'sub'
        elif 'delins' in change:
            operation = 'delins'
        elif 'insdel' in change:
            operation = 'insdel'
        elif 'ins' in change:
            operation = 'ins'
        elif 'dup' in change:
            operation = 'dup'
        else:
            operation = 'unknown'

        return Variation(
            acc_id=acc_id,
            start=int(start.base),
            end=int(end.base),
            ref=ref,
            alt=alt,
            operation=operation,
            source=source,
            coordinate_type=hgvs_variation.type
        )

class VariationLocation(models.Model):
    seq_id    = models.CharField(max_length=30)
    start     = models.IntegerField(null=True)
    end       = models.IntegerField(null=True)
    source    = models.CharField(max_length=30, default='Unknown')
    hgvs      = models.CharField(max_length=200, null=True)

class Pathology(models.Model):
    name = models.CharField(max_length=200)
    # Path to pickled model
    prediction_model = models.CharField(max_length=200)
    variations = models.ManyToManyField(VariationLocation)
    # Precision reported by classification_report function
    model_precision = models.FloatField()